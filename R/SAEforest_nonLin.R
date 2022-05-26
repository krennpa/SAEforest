#' Main function for the estimation of domain-level nonlinear indicators with MERFs under unit-level data
#'
#' This function enables the use of Mixed Effects Random Forests (MERFs)
#' for applications of Small Area Estimation (SAE). Unit-level survey-sample and unit-level
#' covariate data on predictive covariates is required to produce reliable estimates of various
#' disaggregated economic and inequality indicators. The set of predefined indicators includes
#' the mean, median, quantiles (10\%, 25\%, 75\% and 90\%), the head count ratio, the poverty gap,
#' the Gini coefficient and the quintile share ratio. The MERF algorithm is an algorithmic procedure
#' reminiscent of an EM-algorithm (see Details). Overall, the function serves as a coherent framework
#' for the estimation of point estimates and if requested for assessing the uncertainty of the
#' estimates. Methodological details are provided by Krennmair & Schmid (202X) and following
#' examples showcase potential applications.
#'
#' @param Y Continuous input value of target variable.
#' @param X Matrix or data.frame of predictive covariates.
#' @param dName Character specifying the name of domain identifier, for which random intercepts
#' are modeled.
#' @param smp_data data.frame of survey sample data including the specified elements of \code{Y} and
#' \code{X}.
#' @param pop_data data.frame of unit-level population or census level covariate data for
#' covariates \code{X}. Please note that the column names of predictive covariates must match
#' column names of \code{smp_data}. This holds especially for the name of the domain identifier.
#' @param mse Character input specifying the type of uncertainty estimates. Available options are:
#' (i) "none" if only point estimates are requested,
#' (ii) "nonparametric" following the mse bootstrap procedure proposed by Krennmair & Schmid (202X) or
#' (iii) "wild" following the mse bootstrap procedure proposed by Krennmair & Schmid (202X). Defaults to "none".
#' @param importance Variable importance mode processed by the
#' random forest from \pkg{ranger}. Must be 'none', 'impurity', 'impurity_corrected',
#' 'permutation'. Defaults to "none". If you wish to produce informative plots with the generic function
#' \code{\link{plot}}, set \code{importance} not to 'none'. For further details see \link[ranger]{ranger}.
#' @param initialRandomEffects Numeric value or vector of initial estimate of random effects.
#' Defaults to 0.
#' @param ErrorTolerance Numeric value to monitor the MERF algorithm's convergence. Defaults to 1e-04.
#' @param MaxIterations Numeric value specifying the maximal amount of iterations for the
#' MERF algorithm. Defaults to 25.
#' @param B Number of bootstrap replications for mse estimation procedure proposed by
#' Krennmair et al. (202X). Defaults to 100.
#' @param B_adj Number of bootstrap replications for the adjustment of residual variance proposed
#' by Mendez and Lohr (2001). Defaults to 100.
#' @param na.rm Logical. Whether missing values should be removed. Defaults to \code{TRUE}.
#' @param ... Additional parameters are directly passed to the random forest \link[ranger]{ranger}.
#' Most important parameters are for instance \code{mtry} (number of variables to possibly split at
#' in each node), or \code{num.trees} (number of trees). For further details on possible parameters
#' see \link[ranger]{ranger} and the example below.
#' @param threshold Set a custom threshold for indicators, such as the head count ratio. The threshold
#' can be a known numeric value or function of \code{Y}. If the threshold is \code{NULL}, 60 \% of the
#' median of \code{Y} are taken as threshold. Defaults to NULL.
#' @param custom_indicator A list of additional functions containing the indicators to be
#' calculated. These functions must only depend on the target variable \code{Y} and optionally the
#' \code{threshold}. Defaults to \code{NULL}
#' @param smearing Logical input indicating whether a smearing based approach or a MC-based version for
#' point estimates should be obtained. MC should be used if computational constraints do not allow for a
#' smearing based approach. For theoretical details see (WP). Defaults to \code{TRUE}.
#' @param B_MC Number of bootstrap populations to be generated for the MC version for estimating
#' point estimates. Defaults to 100.
#'
#' @return An object of class \code{SAEforest} includes point estimates for disaggregated indicators
#' as well as information on the MERF-model. Optionally corresponding MSE estimates are returned.
#' Several generic functions have methods for the returned object of class \code{SAEforest}. Additionally,
#' the included \code{MERFmodel} object allows the use of generic functions for classes \code{ranger} and
#' \code{\link[lme4]{merMod}}. For a full list and explanation of components and possibilities for objects of class \code{SAEforest},
#' see \code{\link{SAEforestObject}}.
#'
#' @details
#' Mixed effects random forests combine advantages of regression forests (such as implicit model-selection
#' and robustness properties) with the ability to model hierarchical dependencies.
#'
#' The MERF algorithm iteratively optimizes two separate steps: a) the random forest function, assuming the
#' random effects term to be correct and b) estimates the random effects part, assuming the OOB-predictions
#' from the forest to be correct. Overall convergence of the algorithm is monitored by log-likelihood of a
#' joint model of both components. For further details see Krennmair and Schmid or Hajem et. al. (2014).
#'
#' For the estimation of (nonlinear) poverty estimators and/or quantiles, we need information on the
#' area-specific CDF of \code{Y}. We follow a smearing approach originating from Duan (1983) and analyzed within
#' a general unit-level framework for the estimation of SAE means and quantiles by Tzavidis et al. (2010).
#' For further details please see Krennmair et al. (202X). Alternatively to the smearing approach, it is
#' possible to simulate population values of \code{Y} using Monte-Carlo methods. This option is discussed as well
#' and should be used for cases where smearing is not possible due to computational constraints.
#'
#' For the estimation of the MSE, the bootstrap population is built based on a bias-corrected residual
#' variance as proposed Krennmair and Schmid (202X). The bootstrap bias correction follows Mendez and Lohr (2011).
#'
#' Note that the \code{MERFmodel} object is a composition of elements from a random forest of class
#' \code{ranger} and a random effects model of class \code{\link[lme4]{merMod}}. Thus, all generic functions are
#' applicable to corresponding objects. For further details on generic functions see \code{\link[ranger]{ranger}}
#' and \code{\link[lme4]{lmer}} as well as the examples below.
#'
#' @references
#' Krennmair, P. and Schmid, T. (202X). WP 1
#'
#' Krennmair, P. et al. (202X). WP 2
#'
#' Mendez, G. and Lohr, S. (2011) Paper
#'
#' @seealso \code{\link{SAEforestObject}}, \code{\link[ranger]{ranger}}, \code{\link[lme4]{lmer}}
#'
#' @examples
#' \dontrun{
#' #Loading data
#' data("eusilcA_pop")
#' data("eusilcA_smp")
#'
#' income <- eusilcA_smp$eqIncome
#' X_covar <- eusilcA_smp[,-c(1,16,17,18)]
#'
#' #Example 1:
#' #Calculating point estimates and discussing basic generic functions
#'
#' model1 <- SAEforest_nonLin(Y = income, X = X_covar, dName = "district", smp_data = eusilcA_smp,
#'                            pop_data = eusilcA_pop)
#'
#'#example of SAEforest generic
#'summary(model1)
#'
#'#Example 2:
#'#Calculating point estimates and discussing basic generic functions
#'
#'model2 <- SAEforest_nonLin(Y = income, X = X_covar, dName = "district", smp_data = eusilcA_smp,
#'                            pop_data = eusilcA_pop, smearing = FALSE)
#'
#'#example of SAEforest generic
#'summary(model2)
#'
#'
#'#Example 3:
#'#Calculating point + MSE estimates and passing arguments to the forest.
#'#Two additional indicators and function as threshold are added.
#'#Note that B is unrealistically low to improve example speed.
#'
#'model3 <- SAEforest_nonLin(Y = income, X = X_covar, dName = "district", smp_data = eusilcA_smp,
#'                           pop_data = eusilcA_pop, mse = "nonparametric", B = 5, mtry=5,
#'                           num.trees = 100, threshold = function(Y){0.5 * median(Y)},
#'                           custom_indicator = list(my_max = function(Y, threshold){max(Y)},
#'                           my_quant = function(Y, threshold){quantile(Y, probs=c(0.05,0.95))}),
#'                           smearing = FALSE)
#'
#'#example of SAEforest generic:
#'summary(model3)
#'summarize_indicators(model3, MSE = FALSE, CV =TRUE, indicator = c("Gini", "my_max", "my_quant.5%",
#'"my_quant.95%"))
#'}
#'
#' @export

SAEforest_nonLin <- function(Y, X, dName, smp_data, pop_data, smearing = TRUE, mse = "none", importance ="none",
                           initialRandomEffects = 0, ErrorTolerance = 0.0001,
                           MaxIterations = 25, B=100, B_adj=100, B_MC=100, threshold = NULL,
                           custom_indicator =NULL, na.rm = TRUE, ...){

  # ERROR CHECKS OF INPUTS
  #________________________________________

  if(na.rm == TRUE){
    comp_smp <- complete.cases(smp_data)
    smp_data <- smp_data[comp_smp,]
    Y <- Y[comp_smp]
    X <- X[comp_smp,]
    pop_data <- pop_data[complete.cases(pop_data),]
  }

  input_checks_nonLin(Y = Y, X = X, dName = dName, smp_data = smp_data, pop_data =pop_data, smearing =smearing,
                      initialRandomEffects =initialRandomEffects, ErrorTolerance =ErrorTolerance,
                      MaxIterations =MaxIterations, mse=mse, B=B, B_adj=B_adj,B_MC=B_MC, threshold = threshold,
                      importance = importance, custom_indicator = custom_indicator, na.rm =na.rm)

  out_call <- match.call()

  # Make domain variable to character and sort data-sets
  smp_data[[dName]] <- factor(smp_data[[dName]], levels=unique(smp_data[[dName]]))
  pop_data[[dName]] <- factor(pop_data[[dName]], levels=unique(pop_data[[dName]]))

  # Order Data according to factors to ease MSE estimation
  order_in <- order(smp_data[[dName]])
  smp_data <- smp_data[order_in,]
  X <- X[order_in,]
  Y <- Y[order_in]
  pop_data <- pop_data[order(pop_data[[dName]]),]


  # SMEARING
if(smearing == TRUE){
  # Point Estimation
  #________________________________________

  nonLin_preds <- point_nonLin(Y = Y, X = X, dName = dName, threshold = threshold, smp_data = smp_data,
                               pop_data = pop_data, initialRandomEffects = initialRandomEffects,
                               ErrorTolerance = ErrorTolerance, MaxIterations = MaxIterations, importance = importance,
                               custom_indicator=custom_indicator, ...)

  data_specs <- sae_specs(dName = dName,cns = pop_data,smp = smp_data)

  if(mse == "none"){
    result <- list(
      MERFmodel = c(nonLin_preds[[2]], call = out_call, data_specs = list(data_specs), data=list(smp_data)),
      Indicators = sortAlpha(nonLin_preds[[1]],dName =dName),
      MSE_Estimates = NULL,
      AdjustedSD = NULL)

    class(result) <- c("SAEforest_nonLin", "SAEforest")
    return(result)
  }

  # MSE Estimation
  #________________________________________

  if(mse != "none"){
    print(paste("Error SD Bootstrap started:"))
    adj_SD <- adjust_ErrorSD(Y=Y, X=X, smp_data = smp_data, mod = nonLin_preds[[2]], B = B_adj, ...)
    print(paste("MSE Bootstrap with", B,"rounds started:"))
  }

  if(mse == "wild"){
    mse_estims <- MSE_SAEforest_nonLin_wild(Y=Y, X = X, mod=nonLin_preds[[2]], smp_data = smp_data,
                                            pop_data = pop_data, dName = dName, ADJsd = adj_SD , B=B,
                                            threshold = threshold, initialRandomEffects = initialRandomEffects,
                                            ErrorTolerance = ErrorTolerance, MaxIterations = MaxIterations,
                                            custom_indicator =custom_indicator, ...)

    result <- list(
      MERFmodel = c(nonLin_preds[[2]], call = out_call, data_specs = list(data_specs), data=list(smp_data)),
      Indicators = sortAlpha(nonLin_preds[[1]],dName =dName),
      MSE_Estimates = sortAlpha(mse_estims,dName =dName),
      AdjustedSD = adj_SD)

    class(result) <- c("SAEforest_nonLin", "SAEforest")
    return(result)
  }

  if(mse == "nonparametric"){
    mse_estims <- MSE_SAEforest_nonLin_REB(Y=Y, X = X, mod=nonLin_preds[[2]], smp_data = smp_data,
                                           pop_data = pop_data, dName = dName, ADJsd = adj_SD , B=B,
                                           threshold = threshold, initialRandomEffects = initialRandomEffects,
                                           ErrorTolerance = ErrorTolerance, MaxIterations = MaxIterations,
                                           custom_indicator =custom_indicator, ...)

    result <- list(
      MERFmodel = c(nonLin_preds[[2]], call = out_call, data_specs = list(data_specs), data=list(smp_data)),
      Indicators = sortAlpha(nonLin_preds[[1]],dName =dName),
      MSE_Estimates = sortAlpha(mse_estims,dName =dName),
      AdjustedSD = adj_SD)

    class(result) <- c("SAEforest_nonLin", "SAEforest")
    return(result)
  }

}


if(smearing == FALSE){

  nonLin_preds <- point_MC_nonLin(Y = Y, X = X, dName = dName, threshold = threshold, smp_data = smp_data,
                                  pop_data = pop_data, initialRandomEffects = initialRandomEffects,
                                  ErrorTolerance = ErrorTolerance, MaxIterations = MaxIterations, importance = importance,
                                  custom_indicator=custom_indicator, B_point = B_MC, ...)

  data_specs <- sae_specs(dName = dName,cns = pop_data,smp = smp_data)

  if(mse == "none"){
    result <- list(
    MERFmodel = c(nonLin_preds[[2]], call = out_call, data_specs = list(data_specs), data=list(smp_data)),
    Indicators = sortAlpha(nonLin_preds[[1]],dName =dName),
    MSE_Estimates = NULL,
    AdjustedSD = NULL)

    class(result) <- c("MC_MERF_nonLin", "SAEforest")
    return(result)
  }

# MSE Estimation
#________________________________________

if(mse != "none"){
  print(paste("Error SD Bootstrap started:"))
  adj_SD <- adjust_ErrorSD(Y=Y, X=X, smp_data = smp_data, mod = nonLin_preds[[2]], B = B_adj, ...)
  print(paste("MSE Bootstrap with", B,"rounds started:"))
}

if(mse == "wild"){
  mse_estims <- MSE_MC_MERF_nonLin_wild(Y=Y, X = X, mod=nonLin_preds[[2]], smp_data = smp_data,
                                        pop_data = pop_data, dName = dName, ADJsd = adj_SD , B=B, B_point =B_MC,
                                        threshold = threshold, initialRandomEffects = initialRandomEffects,
                                        ErrorTolerance = ErrorTolerance, MaxIterations = MaxIterations,
                                        custom_indicator =custom_indicator, ...)

  result <- list(
    MERFmodel = c(nonLin_preds[[2]], call = out_call, data_specs = list(data_specs), data=list(smp_data)),
    Indicators = sortAlpha(nonLin_preds[[1]],dName =dName),
    MSE_Estimates = sortAlpha(mse_estims,dName =dName),
    AdjustedSD = adj_SD)

  class(result) <- c("MC_MERF_nonLin", "SAEforest")
  return(result)
}

if(mse == "nonparametric"){
  mse_estims <- MSE_MC_MERF_nonLin_REB(Y=Y, X = X, mod=nonLin_preds[[2]], smp_data = smp_data,
                                       pop_data = pop_data, dName = dName, ADJsd = adj_SD , B=B, B_point =B_MC,
                                       threshold = threshold, initialRandomEffects = initialRandomEffects,
                                       ErrorTolerance = ErrorTolerance, MaxIterations = MaxIterations,
                                       custom_indicator =custom_indicator, ...)

  result <- list(
    MERFmodel = c(nonLin_preds[[2]], call = out_call, data_specs = list(data_specs), data=list(smp_data)),
    Indicators = sortAlpha(nonLin_preds[[1]],dName =dName),
    MSE_Estimates = sortAlpha(mse_estims,dName =dName),
    AdjustedSD = adj_SD)

  class(result) <- c("MC_MERF_nonLin", "SAEforest")
  return(result)
}
}
}










