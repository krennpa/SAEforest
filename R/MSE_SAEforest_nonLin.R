# MSE bootstrap for nonlinear indicators under unit-level covariate data ------------------

MSE_SAEforest_nonLin <- function(Y,
                                 X,
                                 dName,
                                 threshold,
                                 smp_data,
                                 mod,
                                 ADJsd,
                                 pop_data,
                                 B = 100,
                                 B_point,
                                 initialRandomEffects = 0,
                                 ErrorTolerance = 0.0001,
                                 MaxIterations = 25,
                                 custom_indicator,
                                 wild,
                                 MC,
                                 aggregate_to,
                                 ...) {

  rand_struc <- paste0(paste0("(1|", dName), ")")
  domains <- t(unique(pop_data[dName]))
  in_samp <- domains %in% t(unique(smp_data[dName]))
  N_i <- as.numeric(table(pop_data[[dName]]))

  pred_vals <- predict(mod$Forest, pop_data)$predictions
  pred_mat <- matrix(pred_vals, nrow = length(pred_vals), ncol = B)

  # Prepare data for sampling
  if (wild == TRUE) {
    ran_obj <- ran_comp_wild(Y = Y, smp_data = smp_data, mod = mod, ADJsd = ADJsd, dName = dName)
    ran_effs <- ran_obj$ran_effs
    forest_res <- ran_obj$forest_res

    sample_e <- function(x) {
      wild_errors(x = x, mod = mod, smp_data = smp_data, forest_res = forest_res)
    }
  }

  if (wild == FALSE) {
    ran_obj <- ran_comp(Y = Y, smp_data = smp_data, mod = mod, ADJsd = ADJsd, dName = dName)
    ran_effs <- ran_obj$ran_effs
    forest_res <- ran_obj$forest_res
    smp_data <- ran_obj$smp_data

    sample_e <- function(x) {
      sample(forest_res, size = sum(N_i), replace = TRUE)
    }
  }

  sample_ui <- function(x) {
    rep(sample(ran_effs,
      size = length(N_i),
      replace = TRUE
    ), N_i)
  }

  u_i <- apply(pred_mat, 2, sample_ui)

  smp_data$forest_res <- NULL

  # combine to y_star
  mu_ij <- pred_mat + u_i
  e_ij <- matrix(NA, nrow = length(pred_vals), ncol = B)
  e_ij <- apply(mu_ij, 2, sample_e)

  y_star <- mu_ij + e_ij

  # get tau_star
  y_star_L <- split(y_star, col(y_star))
  thresh_L <- sapply(y_star_L, function(x) {
    get_thresh(x, threshold = threshold)
  }, simplify = FALSE)
  y_star_L <- Map(cbind, "y_star" = y_star_L, "thresh" = thresh_L)

  if(is.null(aggregate_to)){
    my_agg <- function(x) {
    tapply(x[, 1], pop_data[[dName]], calc_indicat, threshold = unique(x[, 2]), custom = custom_indicator)
    }
  } else{
    my_agg <- function(x) {
      tapply(x[, 1], pop_data[[aggregate_to]], calc_indicat, threshold = unique(x[, 2]), custom = custom_indicator)
    }
  }

  tau_star <- sapply(y_star_L, my_agg, simplify = FALSE)

  if(is.null(aggregate_to)){
    comb <- function(x) {
    matrix(unlist(x), nrow = length(N_i), byrow = TRUE)
    }
  } else{
    N_i_agg <- as.numeric(table(pop_data[[aggregate_to]]))
    comb <- function(x) {
      matrix(unlist(x), nrow = length(N_i_agg), byrow = TRUE)
    }
  }

  tau_star <- sapply(tau_star, comb, simplify = FALSE)

  # get bootstrap samples
  boots_sample <- vector(mode = "list", length = B)

  for (i in 1:B) {
    pop_data$y_star <- y_star[, i]
    boots_sample[[i]] <- sample_select(pop_data, smp = smp_data, dName = dName)
  }

  # uses sample to estimate tau_b
  if (MC == TRUE) {
    my_estim_f <- function(x) {
      point_MC_nonLin(
        Y = x$y_star, X = x[, colnames(X)], dName = dName, threshold = threshold, smp_data = x, pop_data = pop_data,
        initialRandomEffects = initialRandomEffects, ErrorTolerance = ErrorTolerance, B_point = B_point,
        MaxIterations = MaxIterations, custom_indicator = custom_indicator, aggregate_to = aggregate_to, ...
      )[[1]][, -1]
    }
  }

  if (MC == FALSE) {
    my_estim_f <- function(x) {
      point_nonLin(
        Y = x$y_star, X = x[, colnames(X)], dName = dName, threshold = threshold, smp_data = x, pop_data = pop_data,
        initialRandomEffects = initialRandomEffects, ErrorTolerance = ErrorTolerance,
        MaxIterations = MaxIterations, custom_indicator = custom_indicator, aggregate_to = aggregate_to, ...
      )[[1]][, -1]
    }
  }

  tau_b <- pbapply::pbsapply(boots_sample, my_estim_f, simplify = FALSE)

  mean_square <- function(x, y) {
    (x - y)^2
  }

  Mean_square_B <- mapply(mean_square, tau_b, tau_star, SIMPLIFY = FALSE)

  MSE_estimates <- Reduce("+", Mean_square_B) / length(Mean_square_B)

  if(is.null(aggregate_to)){
    MSE_estimates_out <- data.frame(unique(pop_data[dName]), MSE_estimates)
    colnames(MSE_estimates_out) <- c(dName, colnames(MSE_estimates))
  } else{
    MSE_estimates_out <- data.frame(unique(pop_data[aggregate_to]), MSE_estimates)
    colnames(MSE_estimates_out) <- c(aggregate_to, colnames(MSE_estimates))
  }

  rownames(MSE_estimates_out) <- NULL

  # __________________________
  return(MSE_estimates_out)
}
