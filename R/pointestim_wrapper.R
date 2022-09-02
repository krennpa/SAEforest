# Wrapper Functions for Point Estimates
#
# Point estimates for mean and unit-level data---------------------------------------------

point_mean <- function(Y,
                       X,
                       dName,
                       smp_data,
                       pop_data,
                       initialRandomEffects,
                       ErrorTolerance,
                       MaxIterations,
                       importance = "none",
                       ...) {

  random <- paste0(paste0("(1|", dName), ")")

  unit_model <- MERFranger(
    Y = Y,
    X = X,
    random = random,
    data = smp_data,
    initialRandomEffects = initialRandomEffects,
    ErrorTolerance = ErrorTolerance,
    MaxIterations = MaxIterations,
    importance = importance,
    ...
  )

  unit_preds <- predict(unit_model$Forest, pop_data)$predictions +
    predict(unit_model$EffectModel, pop_data, allow.new.levels = TRUE)

  unit_preds_ID <- cbind(pop_data[dName], unit_preds)

  f0 <- as.formula(paste0("unit_preds ", " ~ ", dName))

  mean_preds <- aggregate(f0,
    data = unit_preds_ID,
    FUN = mean
  )
  colnames(mean_preds) <- c(dName, "Mean")

  out_ob <- vector(mode = "list", length = 2)

  out_ob[[1]] <- mean_preds
  out_ob[[2]] <- unit_model

  return(out_ob)
}


# Point estimates for nonlinear indicators and unit-level data (smearing)------------------

point_nonLin <- function(Y,
                         X,
                         dName,
                         threshold,
                         smp_data,
                         pop_data,
                         initialRandomEffects,
                         ErrorTolerance,
                         MaxIterations,
                         importance = "none",
                         custom_indicator,
                         ...) {

  random <- paste0(paste0("(1|", dName), ")")
  domains <- names(table(pop_data[[dName]]))
  popSize <- as.numeric(table(pop_data[[dName]]))

  thresh <- get_thresh(Y, threshold = threshold)

  unit_model <- MERFranger(
    Y = Y,
    X = X,
    random = random,
    data = smp_data,
    initialRandomEffects = initialRandomEffects,
    ErrorTolerance = ErrorTolerance,
    MaxIterations = MaxIterations,
    importance = importance,
    ...
  )

  unit_preds <- predict(unit_model$Forest, pop_data)$predictions +
    predict(unit_model$EffectModel, pop_data, allow.new.levels = TRUE)


  #  smearing step
  smear_list <- vector(mode = "list", length = length(domains))

  for (i in seq_along(domains)) {
    smear_i <- matrix(rep(unit_model$OOBresiduals, popSize[i]), nrow = popSize[i], ncol = length(unit_model$OOBresiduals), byrow = TRUE)
    smear_i <- smear_i + unit_preds[pop_data[[dName]] == domains[i]]
    val_i <- c(smear_i)
    smear_list[[i]] <- calc_indicat(val_i, threshold = thresh, custom = custom_indicator)
  }

  indicators <- do.call(rbind.data.frame, smear_list)
  indicators_out <- cbind(domains, indicators)
  names(indicators_out)[1] <- dName

  out_ob <- vector(mode = "list", length = 2)

  out_ob[[1]] <- indicators_out
  out_ob[[2]] <- unit_model

  return(out_ob)
}


# Point estimates for mean and aggregated covariate data-----------------------------------

point_meanAGG <- function(Y,
                          X,
                          dName,
                          smp_data,
                          Xpop_agg,
                          initialRandomEffects,
                          ErrorTolerance,
                          MaxIterations,
                          OOsample_obs,
                          ADDsamp_obs,
                          w_min,
                          wSet = NULL,
                          importance = "none",
                          verbose = TRUE,
                          samp_seed = 0712,
                          ...) {

  random <- paste0(paste0("(1|", dName), ")")
  groupNames <- as.character(unique(smp_data[[dName]]))
  groupNamesCens <- as.character(unique(Xpop_agg[[dName]]))
  OOsamp <- groupNamesCens[!groupNamesCens %in% groupNames]


  # similarity for out-of-sample
  similarXcens <- Xpop_agg[, colnames(Xpop_agg) != dName]
  rownames(similarXcens) <- groupNamesCens
  simXcensMatrix <- as.matrix(dist(similarXcens))
  diag(simXcensMatrix) <- NA

  unit_model <- MERFranger(
    Y = Y,
    X = X,
    random = random,
    data = smp_data,
    initialRandomEffects = initialRandomEffects,
    ErrorTolerance = ErrorTolerance,
    MaxIterations = MaxIterations,
    importance = importance,
    ...
  )

  unit_preds <- predict(unit_model$Forest, smp_data)$predictions
  u_ij <- predict(unit_model$EffectModel, smp_data, allow.new.levels = TRUE)

  smp_data$forest_preds <- unit_preds + u_ij
  smp_data$u_ij <- u_ij

  joint_smp_data <- smp_data

  # order vars by simularity
  if (is.null(wSet)) {
    wSet <- names(sort(ranger::importance(unit_model$Forest), decreasing = TRUE))
  }

  wSet <- wSet[wSet %in% names(Xpop_agg)]

  if (length(OOsamp) != 0) {
    sim_groups <- apply(simXcensMatrix[!(groupNamesCens %in% OOsamp), OOsamp], 2, FUN = which.min)
    sim_groups <- groupNamesCens[!(groupNamesCens %in% OOsamp)][sim_groups]

    # sampling and incorporating OOsample data
    smp_oos <- vector(mode = "list", length = length(OOsamp))

    for (i in seq(length(OOsamp))) {
      samp_from <- smp_data[as.character(smp_data[[dName]]) == sim_groups[i], ]
      set.seed(samp_seed)
      return_oos <- dplyr::sample_n(samp_from, OOsample_obs, replace = TRUE)
      return_oos[dName] <- OOsamp[i]
      smp_oos[[i]] <- return_oos
    }

    OOs_smp_data <- do.call(rbind.data.frame, smp_oos)

    # out-of-sample observations
    unit_preds_add <- predict(unit_model$Forest, OOs_smp_data)$predictions
    u_ij <- 0

    OOs_smp_data$forest_preds <- unit_preds_add
    OOs_smp_data$u_ij <- u_ij

    joint_smp_data <- rbind(smp_data, OOs_smp_data)
  }

  # find weights and adjust for failure
  smp_weightsIncluded <- vector(mode = "list", length = length(groupNamesCens))
  smp_weightsNames <- vector(mode = "list", length = length(groupNamesCens))

  for (i in groupNamesCens) {
    pos <- which(i == groupNamesCens)

    X_input_elm <- as.matrix(joint_smp_data[wSet])[joint_smp_data[dName] == i, ]
    mu_input_elm <- as.matrix(Xpop_agg[wSet])[Xpop_agg[dName] == i, ]

    X_input_elm <- X_input_elm[, colSums(X_input_elm != 0) > 0]
    mu_input_elm <- mu_input_elm[colnames(X_input_elm)]

    ELMweight <- elm_wrapper(X_input_elm, mu_input_elm)
    sum_w <- round(sum(ELMweight$prob), digits = 7)

    w_smp_data <- joint_smp_data[joint_smp_data[dName] == i, ]

    if (sum_w == 1) {
      w_smp_data$weights <- ELMweight$prob
      rownames(w_smp_data) <- NULL
      smp_weightsIncluded[[pos]] <- w_smp_data
      smp_weightsNames[[pos]] <- colnames(X_input_elm)
    } else {
      mod_smp_data <- joint_smp_data[joint_smp_data[dName] == i, ]
      rownames(mod_smp_data) <- NULL

      if (ADDsamp_obs != 0) {
        similarXcens <- Xpop_agg[, colnames(Xpop_agg) != dName]
        rownames(similarXcens) <- groupNamesCens
        simXcensMatrix <- as.matrix(dist(similarXcens))
        diag(simXcensMatrix) <- NA

        sim_group <- which.min(simXcensMatrix[, pos])
        samp_add <- joint_smp_data[joint_smp_data[dName] == groupNamesCens[sim_group], ]
        set.seed(samp_seed)
        return_add <- dplyr::sample_n(samp_add, ADDsamp_obs, replace = TRUE)
        return_add[dName] <- i

        mod_smp_data <- rbind(
          joint_smp_data[joint_smp_data[dName] == i, ],
          return_add
        )

        rownames(mod_smp_data) <- NULL
      }

      # recalculation
      X_input_elm <- as.matrix(mod_smp_data[wSet])
      mu_input_elm <- as.matrix(Xpop_agg[wSet])[Xpop_agg[dName] == i, ]

      X_input_elm <- X_input_elm[, colSums(X_input_elm != 0) > 0]
      mu_input_elm <- mu_input_elm[colnames(X_input_elm)]

      ELMweight <- elm_wrapper(X_input_elm, mu_input_elm)
      sum_w <- round(sum(ELMweight$prob), digits = 7)

      if (sum_w == 1) {
        mod_smp_data$weights <- ELMweight$prob
        rownames(mod_smp_data) <- NULL
        smp_weightsIncluded[[pos]] <- mod_smp_data
        smp_weightsNames[[pos]] <- colnames(X_input_elm)
      } else {
        sum_w <- 0

        while ((sum_w != 1) & (dim(X_input_elm)[2] > w_min)) {
          X_input_elm <- X_input_elm[, -dim(X_input_elm)[2]]
          mu_input_elm <- mu_input_elm[-dim(X_input_elm)[2]]

          ELMweight <- elm_wrapper(X_input_elm, mu_input_elm)
          sum_w <- round(sum(ELMweight$prob), digits = 7)
        }
      }

      if (sum_w == 1) {
        mod_smp_data$weights <- ELMweight$prob
        rownames(mod_smp_data) <- NULL
        smp_weightsIncluded[[pos]] <- mod_smp_data
        smp_weightsNames[[pos]] <- colnames(X_input_elm)
      } else {
        if (verbose == TRUE) {
          message(paste("Calculation of weights failed for area:", i))
        }
        rownames(w_smp_data) <- NULL
        w_smp_data$weights <- 1 / dim(w_smp_data)[1]
        smp_weightsIncluded[[pos]] <- w_smp_data
        smp_weightsNames[[pos]] <- NA
      }
    }
  }

  final_smp_data <- do.call(dplyr::bind_rows, smp_weightsIncluded)

  # estimate final means
  final_smp_data$W_mean <- with(final_smp_data, forest_preds * weights)
  f0 <- as.formula(paste("W_mean", paste(dName), sep = " ~ "))
  Mean_preds <- aggregate(f0, data = final_smp_data, FUN = sum)
  colnames(Mean_preds)[2] <- c("Mean")

  # Prepare return object
  return(list(
    Indicators = Mean_preds,
    MERFmodel = unit_model,
    ModifiedSet = final_smp_data,
    ADDsamp_obs = ADDsamp_obs,
    OOsample_obs = OOsample_obs,
    wSet = wSet, w_min = w_min,
    wAreaInfo = smp_weightsNames
  ))
}


# Point estimates for nonlinear indicators and unit-level data (MC version) ---------------

point_MC_nonLin <- function(Y,
                            X,
                            dName,
                            threshold,
                            smp_data,
                            pop_data,
                            initialRandomEffects,
                            B_point,
                            ErrorTolerance,
                            MaxIterations,
                            importance = "none",
                            custom_indicator,
                            ...) {

  domains <- names(table(pop_data[[dName]]))
  random <- paste0(paste0("(1|", dName), ")")

  popSize_N <- data.frame(table(pop_data[[dName]]))
  popSize_n <- data.frame(table(smp_data[[dName]]))
  colnames(popSize_N) <- c(dName, "N_i")
  colnames(popSize_n) <- c(dName, "n_i")
  popSize <- dplyr::left_join(popSize_N, popSize_n, by = dName)
  popSize[, 3][is.na(popSize[, 3])] <- 0

  thresh <- get_thresh(Y, threshold = threshold)

  unit_model <- MERFranger(
    Y = Y,
    X = X,
    random = random,
    data = smp_data,
    initialRandomEffects = initialRandomEffects,
    ErrorTolerance = ErrorTolerance,
    MaxIterations = MaxIterations,
    importance = importance,
    ...
  )

  unit_preds <- predict(unit_model$Forest, pop_data)$predictions +
    predict(unit_model$EffectModel, pop_data, allow.new.levels = TRUE)

  # preparing data for MC step
  ran_obj <- ran_comp(Y = Y, smp_data = smp_data, mod = unit_model, ADJsd = unit_model$ErrorSD, dName = dName)
  ran_effs <- ran_obj$ran_effs
  forest_res <- ran_obj$forest_res
  smp_data <- ran_obj$smp_data

  pred_val <- matrix(unit_preds,
    ncol = B_point,
    nrow = length(unit_preds), byrow = FALSE
  )

  y_hat <- pred_val + sample(forest_res, size = length(pred_val), replace = TRUE)

  gamm_i <- (unit_model$RanEffSD^2) / (unit_model$RanEffSD^2 + (unit_model$ErrorSD^2 / popSize$n_i))

  v_i <- apply(y_hat, 2, function(x) {
    rep(sample(ran_effs, size = length(popSize$N_i), replace = TRUE), popSize$N_i)
  })
  gamm_1 <- rep((1 - gamm_i), popSize$N_i)
  y_star <- y_hat + gamm_1 * v_i

  indi_agg <- rep(1:length(popSize$N_i), popSize$N_i)
  my_agg <- function(x) {
    tapply(x, indi_agg, calc_indicat, threshold = thresh, custom = custom_indicator)
  }
  tau_star <- apply(y_star, my_agg, MARGIN = 2, simplify = FALSE)

  col_names <- colnames(tau_star[[1]]$`1`)

  comb <- function(x) {
    matrix(unlist(x), nrow = length(popSize$N_i), byrow = TRUE)
  }
  tau_star <- sapply(tau_star, comb, simplify = FALSE)

  indicators <- Reduce("+", tau_star) / length(tau_star)
  colnames(indicators) <- col_names

  # preparing final output
  result <- list(
    Indicators = cbind(popSize[dName], indicators),
    model = unit_model
  )
  return(result)
}
