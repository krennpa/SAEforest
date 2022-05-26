# Main function for means using MERFs with aggregated population data
#
# This function facilitates the use of Mixed Effects Random Forests (MERFs)
# for applications of Small Area Estimation (SAE). Unit-level survey-sample and
# aggregated covariate data on predictive covariates is required to produce reliable
# estimates of domain-specific means. The MERF algorithm is an algorithmic procedure
# reminiscent of an EM-algorithm (see Details). Overall, the function serves as a coherent
# framework for the estimation of point-estimates and if requested for assessing the uncertainty
# of the estimates. Methodological details are provided by Krennmair et al. (202X) and following
# examples showcase potential applications.
#
# @param Y metric input value of target variable.
# @param X matrix or data.frame of predictive covariates.
# @param dName character specifying the name of domain identifier, for which random intercepts
# are modeled.
# @param smp_data data.frame of survey sample data including the specified elements of \code{Y} and
# \code{X}.
# @param Xpop_agg data.frame of aggregated population or census level covariate data for
# covariates \code{X}. Please note that the column names of predictive covariates must match
# column names of \code{smp_data}. This holds especially for the name of the domain identifier.
# @param mse character input specifying the type of uncertainty estimates. Available options are:
# (i) "nonparametric" following the mse boostrap procedure proposed by Krennmair et al. (202X) or
# (ii) "none" if only point estimates are requested. Defaults to "none".
# @param importance variable importance mode processed by the
# random forest from the \pkg{ranger}. Must be one of the following option: (i)'impurity', (ii) 'impurity_corrected'
# or (iii) 'permutation'. Defaults to 'impurity'. variable importance is used to rank covariate information
# in the process of finding suitable calibration weights (see details). For further information regarding
# the measures of importance see \link[ranger]{ranger}.
# @param initialRandomEffects numeric value or vector of initial estimate of random effects.
# Defaults to 0.
# @param ErrorTolerance numeric value to monitor the MERF algorithm's convergence. Defaults to 1e-04.
# @param MaxIterations numeric value specifying the maximal amount of iterations for the
# MERF algorithm. Defaults to 25.
# @param B numeric number of bootstrap replications for mse estimation procedure proposed by
# Krennmair et al. (202X). Defaults to 100.
# @param B_adj numeric number of bootstrap replications for the adjustment of residual variance. Defaults to 100.
# @param na.rm logical. Whether missing values should be removed. Defaults to TRUE.
# @param ... additional parameters are directly passed to the random forest \link[ranger]{ranger}.
# Most important parameters are for instance mtry (number of variables to possibly split at
# in each node), or num.tree (number of trees). For further details on possible parameters
# see \link[ranger]{ranger} and the example below.
# @param popnsize data.frame, which is only needed if an mse is requested. comprising information of
# population size of domains. Please note that the name of the domain identifier must match the column
# name of \code{smp_data}.
# @param OOsample_obs Number of Out-of-sample observations taken from the closest area for potentially unsampled
# areas. Default set to 25.
# @param ADDsamp_obs Number of Out-of-sample observations taken from the closest area if first iteration for the
# calculation of calibration weights fails. Default is set to 0.
# @param w_min Minimal number of covariates from which informative weights are calculated. Defaults to 3.

# @return An object of class "SAEforest" always includes point estimates for disaggregated indicators
# as well as information on the MERF-model. Optionally corresponding MSE estimates are returned.
# Several generic functions have methods for the returned object of class "SAEforest". Additionally,
# the included \code{MERFmodel} object allows the use of generic functions for classes "ranger" and
# "merMod". For a full list and explanation of components and possibilities for objects of class "SAEforest",
# see \code{\link{SAEforestObject}}.
#
# @details Mixed effects random forests combine advantages of regression forests (such as implicit
# model-selection and robustness properties) with the ability to model hierarchical dependencies.
#
# The MERF algorithm iteratively optimizes two separate steps: a) the random forest function, assuming
# the random effects term to be correct and b) estimates the random effects part, assuming the
# OOB-predictions from the forest to be correct. Overall convergence of the algorithm is monitored
# by log-likelihood of a joint model of both components. For further details see Krennmair and Schmid
# or Hajem et. al. (2014).
#
# In contrast to \code{\link{SAEforest_mean}} which requires covariate micro-data, this function adaptively
# incorporates auxiliary information through calibration-weights based on empirical likelihood for the
# estimation of area-level means. See methodlogical details in Krennmair et al. (202X).
#
# For the estimation of the MSE, a bias-corrected residual variance as proposed Krennmair and Schmid
# (202X) is required. The bootstrap bias correction follows Mendez and Lohr (2011).

# Note that the \code{MERFmodel} object is a composition of elements from a random forest of class
# 'ranger' and a random effects model of class 'lmerMod'.  Thus, all generic functions applicable
# to objects of classes 'ranger' and 'lmerMod' can be used on these elements. For furhter details
# on generic functions see \code{\link[ranger]{ranger}} and \code{\link[lme4]{lmer}} as well as
# the examples below.
#
# @seealso \code{\link{SAEforestObject}}, \code{\link[ranger]{ranger}}, \code{\link[lme4]{lmer}}
#
# @references
# Krennmair, P. and Schmid, T. (202X). WP 1
#
# Mendez, G. and Lohr, S. (2011) Paper
#
# @examples
# \dontrun{#Loading data
# data("eusilcA_popAgg")
# data("eusilcA_smp")
# data("popNsize")
#
# income <- eusilcA_smp$eqIncome
# X_covar <- eusilcA_smp[,-c(1,16,17,18)]
#
##Example 1:
##Calculating point-estimates and discussing basic generic functions
#
# model1 <- SAEforest_meanAGG(Y = income, X = X_covar, dName = "district",
#                             smp_data = eusilcA_smp, Xpop_agg = eusilcA_popAgg)
#
##SAEforest generics:
# summary(model1)
#
#
# #Example 2:
# #Calculating point + MSE estimates and passing arguments to the forest
#
# model2 <- SAEforest_meanAGG(Y = income, X = X_covar, dName = "district",
#                             smp_data = eusilcA_smp, Xpop_agg = eusilcA_popAgg,
#                             mse = "nonparametric", popnsize = popNsize,
#                             B = 25, mtry=5, num.trees = 100)
#
##SAEforest generics:
#summary(model2)
#summarize_indicators(model2, MSE = TRUE, CV =TRUE)
#}
#
# @export

SAEforest_meanAGG <- function(Y, X, dName, smp_data, Xpop_agg, mse = "none", importance = "impurity",
                              OOsample_obs = 25, ADDsamp_obs=0, w_min=3, initialRandomEffects = 0,
                              ErrorTolerance = 0.0001, MaxIterations = 25, B=100, B_adj = 100, popnsize,
                              na.rm =TRUE, ...){

  # ERROR CHECKS OF INPUTS
  #________________________________________

  if(na.rm == TRUE){
    comp_smp <- complete.cases(smp_data)
    smp_data <- smp_data[comp_smp,]
    Y <- Y[comp_smp]
    X <- X[comp_smp,]
    Xpop_agg <- Xpop_agg[complete.cases(Xpop_agg),]
  }

  input_checks_meanAGG(Y = Y, X = X, dName = dName, smp_data =smp_data, Xpop_agg =Xpop_agg,
                       initialRandomEffects = initialRandomEffects, ErrorTolerance =ErrorTolerance,
                       MaxIterations =MaxIterations, mse =mse, B = B, popnsize =popnsize,
                       OOsample_obs =OOsample_obs, ADDsamp_obs = ADDsamp_obs, w_min =w_min,
                       importance = importance, na.rm = na.rm)

  out_call <- match.call()


  # Make domain variable to character and sort data-sets
  smp_data[[dName]] <- factor(smp_data[[dName]], levels=unique(smp_data[[dName]]))
  Xpop_agg[[dName]] <- factor(Xpop_agg[[dName]], levels = unique(Xpop_agg[[dName]]))
  if(mse != "none"){
    popnsize[[dName]] <- factor(popnsize[[dName]], levels = unique(Xpop_agg[[dName]]))
  }

  # Order Data according to factors to ease MSE-estimation
  order_in <- order(smp_data[[dName]])
  smp_data <- smp_data[order_in,]
  X <- X[order_in,]
  Y <- Y[order_in]

    # Point Estimation
  #________________________________________
  meanAGG_preds <- point_meanAGG(Y = Y, X = X, dName = dName, smp_data = smp_data, Xpop_agg = Xpop_agg,
                               initialRandomEffects = initialRandomEffects, ErrorTolerance = ErrorTolerance,
                               MaxIterations = MaxIterations, OOsample_obs = OOsample_obs, ADDsamp_obs = ADDsamp_obs, w_min = w_min,
                               importance = importance, ...)

  data_specs <- sae_specs(dName = dName,cns = Xpop_agg, smp = smp_data)


  if(mse == "none"){
    result <- list(
      MERFmodel =  c(meanAGG_preds[[2]], call = out_call, data_specs = list(data_specs), data=list(smp_data)),
      Indicators = meanAGG_preds[[1]],
      MSE_Estimates = NULL,
      AdjustedSD = NULL,
      NrCovar = meanAGG_preds$wAreaInfo)

    class(result) <- c("SAEforest_meanAGG", "SAEforest")
    return(result)
  }

  # MSE Estimation
  #________________________________________

  if(mse != "none"){
    print(paste("Error SD Bootstrap started:"))
    adj_SD <- adjust_ErrorSD(Y=Y, X=X, smp_data = smp_data, mod = meanAGG_preds[[2]], B = B_adj, ...)
    print(paste("Bootstrap with", B,"rounds started"))
  }

  if(mse == "nonparametric"){

    mse_estims <- MSE_SAEforest_aggOOB_wSet(Y=Y, X = X, mod = meanAGG_preds, smp_data = smp_data,
                                            Xpop_agg = Xpop_agg, dName = dName, ADJsd = adj_SD , B=B,
                                            initialRandomEffects = initialRandomEffects,
                                           ErrorTolerance = ErrorTolerance, MaxIterations = MaxIterations,
                                           popnsize = popnsize, ...)

    result <- list(
      MERFmodel =  c(meanAGG_preds[[2]], call = out_call, data_specs = list(data_specs), data=list(smp_data)),
      Indicators = meanAGG_preds[[1]],
      MSE_Estimates = mse_estims,
      AdjustedSD = adj_SD,
      NrCovar = meanAGG_preds$wAreaInfo)

    class(result) <- c("SAEforest_meanAGG", "SAEforest")
    return(result)
  }

}
