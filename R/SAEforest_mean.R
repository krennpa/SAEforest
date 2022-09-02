# Wrapper function for point and MSE estimates for domain-level means

SAEforest_mean <- function(Y,
                           X,
                           dName,
                           smp_data,
                           pop_data,
                           MSE = "none",
                           aggData = FALSE,
                           popnsize = NULL,
                           OOsample_obs = 25,
                           ADDsamp_obs = 0,
                           w_min = 3,
                           importance = "impurity",
                           initialRandomEffects = 0,
                           ErrorTolerance = 0.0001,
                           MaxIterations = 25,
                           B = 100,
                           B_adj = 100,
                           na.rm = TRUE,
                           out_call,
                           ...) {

  # Point and MSE estimates for domain-level means and unit-level data ----------------------

  if (aggData == FALSE) {
    if (na.rm == TRUE) {
      comp_smp <- complete.cases(smp_data)
      smp_data <- smp_data[comp_smp, ]
      Y <- Y[comp_smp]
      X <- X[comp_smp, ]
      pop_data <- pop_data[complete.cases(pop_data[,c(colnames(X),dName)]), ]
    }

    # make domain variable to character and sort data-sets
    smp_data[[dName]] <- factor(smp_data[[dName]], levels = unique(smp_data[[dName]]))
    pop_data[[dName]] <- factor(pop_data[[dName]], levels = unique(pop_data[[dName]]))

    # order data according to factors to ease MSE estimation
    order_in <- order(smp_data[[dName]])
    smp_data <- smp_data[order_in, ]
    X <- X[order_in, ]
    Y <- Y[order_in]
    pop_data <- pop_data[order(pop_data[[dName]]), ]


    # point estimation
    mean_preds <- point_mean(
      Y = Y,
      X = X,
      dName = dName,
      smp_data = smp_data,
      pop_data = pop_data,
      initialRandomEffects = initialRandomEffects,
      ErrorTolerance = ErrorTolerance,
      MaxIterations = MaxIterations,
      importance = importance,
      ...
    )

    data_specs <- sae_specs(dName = dName, cns = pop_data, smp = smp_data)


    if (MSE == "none") {
      result <- list(
        MERFmodel = c(mean_preds[[2]], call = out_call, data_specs = list(data_specs), data = list(smp_data)),
        Indicators = sortAlpha(mean_preds[[1]], dName = dName),
        MSE_Estimates = NULL,
        AdjustedSD = NULL,
        NrCovar = NULL
      )

      class(result) <- c("SAEforest_mean", "SAEforest")
      return(result)
    }

    # MSE estimation

    if (MSE != "none") {
      message(paste("Error SD Bootstrap started:"))
      adj_SD <- adjust_ErrorSD(Y = Y, X = X, smp_data = smp_data, mod = mean_preds[[2]], B = B_adj, ...)
      message(paste("Bootstrap with", B, "rounds started"))
    }

    if (MSE == "nonparametric") {
      MSE_estims <- MSE_SAEforest_mean(
        Y = Y,
        X = X,
        dName = dName,
        smp_data = smp_data,
        mod = mean_preds[[2]],
        ADJsd = adj_SD,
        pop_data = pop_data,
        B = B,
        initialRandomEffects = initialRandomEffects,
        ErrorTolerance = ErrorTolerance,
        MaxIterations = MaxIterations,
        ...
      )

      result <- list(
        MERFmodel = c(mean_preds[[2]], call = out_call, data_specs = list(data_specs), data = list(smp_data)),
        Indicators = sortAlpha(mean_preds[[1]], dName = dName),
        MSE_Estimates = sortAlpha(MSE_estims, dName = dName),
        AdjustedSD = adj_SD,
        NrCovar = NULL
      )

      class(result) <- c("SAEforest_mean", "SAEforest")
      return(result)
    }
  }

  # Point and MSE estimates for domain-level means and aggregated data ----------------------

  if (aggData == TRUE) {
    if (na.rm == TRUE) {
      comp_smp <- complete.cases(smp_data)
      smp_data <- smp_data[comp_smp, ]
      Y <- Y[comp_smp]
      X <- X[comp_smp, ]
      pop_data <- pop_data[complete.cases(pop_data[,c(colnames(X),dName)]), ]
    }

    # make domain variable to character and sort data-sets
    smp_data[[dName]] <- factor(smp_data[[dName]], levels = unique(smp_data[[dName]]))
    pop_data[[dName]] <- factor(pop_data[[dName]], levels = unique(pop_data[[dName]]))
    if (MSE != "none") {
      popnsize[[dName]] <- factor(popnsize[[dName]], levels = unique(pop_data[[dName]]))
    }

    # order data according to factors to ease MSE estimation
    order_in <- order(smp_data[[dName]])
    smp_data <- smp_data[order_in, ]
    X <- X[order_in, ]
    Y <- Y[order_in]

    # point estimation
    mean_preds <- point_meanAGG(
      Y = Y,
      X = X,
      dName = dName,
      smp_data = smp_data,
      Xpop_agg = pop_data,
      initialRandomEffects = initialRandomEffects,
      ErrorTolerance = ErrorTolerance,
      MaxIterations = MaxIterations,
      OOsample_obs = OOsample_obs,
      ADDsamp_obs = ADDsamp_obs,
      w_min = w_min,
      importance = importance,
      ...
    )

    data_specs <- sae_specs(dName = dName, cns = pop_data, smp = smp_data)


    if (MSE == "none") {
      result <- list(
        MERFmodel = c(mean_preds[[2]], call = out_call, data_specs = list(data_specs), data = list(smp_data)),
        Indicators = sortAlpha(mean_preds[[1]], dName = dName),
        MSE_Estimates = NULL,
        AdjustedSD = NULL,
        NrCovar = mean_preds$wAreaInfo
      )

      class(result) <- c("SAEforest_meanAGG", "SAEforest")
      return(result)
    }

    # MSE estimation
    if (MSE != "none") {
      message(paste("Error SD Bootstrap started:"))
      adj_SD <- adjust_ErrorSD(Y = Y, X = X, smp_data = smp_data, mod = mean_preds[[2]], B = B_adj, ...)
      message(paste("Bootstrap with", B, "rounds started"))
    }

    if (MSE == "nonparametric") {
      MSE_estims <- MSE_SAEforest_aggOOB(
        Y = Y,
        X = X,
        mod = mean_preds,
        smp_data = smp_data,
        Xpop_agg = pop_data,
        dName = dName,
        ADJsd = adj_SD,
        B = B,
        initialRandomEffects = initialRandomEffects,
        ErrorTolerance = ErrorTolerance,
        MaxIterations = MaxIterations,
        popnsize = popnsize,
        ...
      )

      result <- list(
        MERFmodel = c(mean_preds[[2]], call = out_call, data_specs = list(data_specs), data = list(smp_data)),
        Indicators = sortAlpha(mean_preds[[1]], dName = dName),
        MSE_Estimates = sortAlpha(MSE_estims, dName = dName),
        AdjustedSD = adj_SD,
        NrCovar = mean_preds$wAreaInfo
      )

      class(result) <- c("SAEforest_meanAGG", "SAEforest")
      return(result)
    }
  }
}
