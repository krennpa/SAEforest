# Wrapper function for point and MSE estimates for domain-level nonlinear indicators

SAEforest_nonLin <- function(Y,
                             X,
                             dName,
                             smp_data,
                             pop_data,
                             smearing = TRUE,
                             MSE = "none",
                             importance = "none",
                             initialRandomEffects = 0,
                             ErrorTolerance = 0.0001,
                             MaxIterations = 25,
                             B = 100,
                             B_adj = 100,
                             B_MC = 100,
                             threshold = NULL,
                             custom_indicator = NULL,
                             na.rm = TRUE,
                             out_call,
                             ...) {

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


  # Point and MSE estimates for domain-level indicators and unit-level data (smearing) ------
  if (smearing == TRUE) {

    # point estimation
    nonLin_preds <- point_nonLin(
      Y = Y,
      X = X,
      dName = dName,
      threshold = threshold,
      smp_data = smp_data,
      pop_data = pop_data,
      initialRandomEffects = initialRandomEffects,
      ErrorTolerance = ErrorTolerance,
      MaxIterations = MaxIterations,
      importance = importance,
      custom_indicator = custom_indicator,
      ...
    )

    data_specs <- sae_specs(dName = dName, cns = pop_data, smp = smp_data)

    if (MSE == "none") {
      result <- list(
        MERFmodel = c(nonLin_preds[[2]], call = out_call, data_specs = list(data_specs), data = list(smp_data)),
        Indicators = sortAlpha(nonLin_preds[[1]], dName = dName),
        MSE_Estimates = NULL,
        AdjustedSD = NULL
      )

      class(result) <- c("SAEforest_nonLin", "SAEforest")
      return(result)
    }

    # MSE estimation
    if (MSE != "none") {
      message(paste("Error SD Bootstrap started:"))
      adj_SD <- adjust_ErrorSD(Y = Y, X = X, smp_data = smp_data, mod = nonLin_preds[[2]], B = B_adj, ...)
      message(paste("MSE Bootstrap with", B, "rounds started:"))
    }

    if (MSE == "wild") {
      MSE_estims <- MSE_SAEforest_nonLin(
        Y = Y,
        X = X,
        mod = nonLin_preds[[2]],
        smp_data = smp_data,
        pop_data = pop_data,
        dName = dName,
        ADJsd = adj_SD,
        B = B,
        threshold = threshold,
        initialRandomEffects = initialRandomEffects,
        ErrorTolerance = ErrorTolerance,
        MaxIterations = MaxIterations,
        custom_indicator = custom_indicator,
        wild = TRUE,
        MC = FALSE,
        B_point = B_MC,
        ...
      )

      result <- list(
        MERFmodel = c(nonLin_preds[[2]], call = out_call, data_specs = list(data_specs), data = list(smp_data)),
        Indicators = sortAlpha(nonLin_preds[[1]], dName = dName),
        MSE_Estimates = sortAlpha(MSE_estims, dName = dName),
        AdjustedSD = adj_SD
      )

      class(result) <- c("SAEforest_nonLin", "SAEforest")
      return(result)
    }

    if (MSE == "nonparametric") {
      MSE_estims <- MSE_SAEforest_nonLin(
        Y = Y,
        X = X,
        mod = nonLin_preds[[2]],
        smp_data = smp_data,
        pop_data = pop_data,
        dName = dName,
        ADJsd = adj_SD,
        B = B,
        threshold = threshold,
        initialRandomEffects = initialRandomEffects,
        ErrorTolerance = ErrorTolerance,
        MaxIterations = MaxIterations,
        custom_indicator = custom_indicator,
        wild = FALSE,
        MC = FALSE,
        B_point = B_MC,
        ...
      )

      result <- list(
        MERFmodel = c(nonLin_preds[[2]], call = out_call, data_specs = list(data_specs), data = list(smp_data)),
        Indicators = sortAlpha(nonLin_preds[[1]], dName = dName),
        MSE_Estimates = sortAlpha(MSE_estims, dName = dName),
        AdjustedSD = adj_SD
      )

      class(result) <- c("SAEforest_nonLin", "SAEforest")
      return(result)
    }
  }


  # Point and MSE estimates for domain-level indicators and unit-level data (MC version) ----

  if (smearing == FALSE) {
    nonLin_preds <- point_MC_nonLin(
      Y = Y,
      X = X,
      dName = dName,
      threshold = threshold,
      smp_data = smp_data,
      pop_data = pop_data,
      initialRandomEffects = initialRandomEffects,
      ErrorTolerance = ErrorTolerance,
      MaxIterations = MaxIterations,
      importance = importance,
      custom_indicator = custom_indicator,
      B_point = B_MC,
      ...
    )

    data_specs <- sae_specs(dName = dName, cns = pop_data, smp = smp_data)

    if (MSE == "none") {
      result <- list(
        MERFmodel = c(nonLin_preds[[2]], call = out_call, data_specs = list(data_specs), data = list(smp_data)),
        Indicators = sortAlpha(nonLin_preds[[1]], dName = dName),
        MSE_Estimates = NULL,
        AdjustedSD = NULL
      )

      class(result) <- c("MC_MERF_nonLin", "SAEforest")
      return(result)
    }

    # MSE estimation
    if (MSE != "none") {
      message(paste("Error SD Bootstrap started:"))
      adj_SD <- adjust_ErrorSD(Y = Y, X = X, smp_data = smp_data, mod = nonLin_preds[[2]], B = B_adj, ...)
      message(paste("MSE Bootstrap with", B, "rounds started:"))
    }

    if (MSE == "wild") {
      MSE_estims <- MSE_SAEforest_nonLin(
        Y = Y,
        X = X,
        mod = nonLin_preds[[2]],
        smp_data = smp_data,
        pop_data = pop_data,
        dName = dName,
        ADJsd = adj_SD,
        B = B,
        threshold = threshold,
        initialRandomEffects = initialRandomEffects,
        ErrorTolerance = ErrorTolerance,
        MaxIterations = MaxIterations,
        custom_indicator = custom_indicator,
        wild = TRUE,
        MC = TRUE,
        B_point = B_MC,
        ...
      )

      result <- list(
        MERFmodel = c(nonLin_preds[[2]], call = out_call, data_specs = list(data_specs), data = list(smp_data)),
        Indicators = sortAlpha(nonLin_preds[[1]], dName = dName),
        MSE_Estimates = sortAlpha(MSE_estims, dName = dName),
        AdjustedSD = adj_SD
      )

      class(result) <- c("MC_MERF_nonLin", "SAEforest")
      return(result)
    }

    if (MSE == "nonparametric") {
      MSE_estims <- MSE_SAEforest_nonLin(
        Y = Y,
        X = X,
        mod = nonLin_preds[[2]],
        smp_data = smp_data,
        pop_data = pop_data,
        dName = dName,
        ADJsd = adj_SD,
        B = B,
        threshold = threshold,
        initialRandomEffects = initialRandomEffects,
        ErrorTolerance = ErrorTolerance,
        MaxIterations = MaxIterations,
        custom_indicator = custom_indicator,
        wild = FALSE,
        MC = TRUE,
        B_point = B_MC,
        ...
      )

      result <- list(
        MERFmodel = c(nonLin_preds[[2]], call = out_call, data_specs = list(data_specs), data = list(smp_data)),
        Indicators = sortAlpha(nonLin_preds[[1]], dName = dName),
        MSE_Estimates = sortAlpha(MSE_estims, dName = dName),
        AdjustedSD = adj_SD
      )

      class(result) <- c("MC_MERF_nonLin", "SAEforest")
      return(result)
    }
  }
}
