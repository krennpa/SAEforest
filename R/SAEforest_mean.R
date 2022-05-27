#' Main function for the estimation of domain-level means with MERFs under unit-level and aggregated data
#'
#' This function enables the use of Mixed Effects Random Forests (MERFs) for applications
#' of Small Area Estimation (SAE). Unit-level sample data and additional unit-level or aggregated
#' data on predictive covariates is required to produce reliable estimates of domain-specific means.
#' The MERF algorithm is an algorithmic procedure reminiscent of an EM-algorithm (see Details).
#' Overall, the function serves as a coherent framework for the estimation of point estimates
#' and if requested for assessing the uncertainty of the estimates. Methodological details are
#' provided by Krennmair & Schmid (202X) and following examples showcase potential applications.
#'
#' @param Y Continuous input value of target variable.
#' @param X Matrix or data.frame of predictive covariates.
#' @param dName Character specifying the name of domain identifier, for which random intercepts
#' are modeled.
#' @param smp_data data.frame of survey sample data including the specified elements of \code{Y} and
#' \code{X}.
#' @param pop_data data.frame of unit-level or aggregated population level data for
#' covariates \code{X}. Please note that the column names of predictive covariates and the domain-level
#' identifier must match the column names of \code{smp_data}. Also note, that if aggregated covariate data
#' is used, the function parameter \code{aggData} must be set to \code{TRUE}.
#' @param MSE Character input specifying the type of uncertainty estimates. Currently available options are:
#' (i) "none" if only point estimates are requested or (ii) "nonparametric" following the MSE boostrap procedure
#' proposed by Krennmair & Schmid (202X) and by Krennmair et al (202X) if \code{aggData = TRUE}.
#' @param importance Variable importance mode processed by the
#' random forest from \pkg{ranger}. Must be one of the following options: (i)'impurity', (ii) 'impurity_corrected'
#' or (iii) 'permutation'. Defaults to 'impurity'. In general, a concept of variable importance is needed
#' for the production of informative plots with the generic function \code{\link{plot}}. In the case of
#' aggregated covariate data, variable importance is needed to rank covariate information in the
#' process of finding suitable calibration weights (see Details). For further information regarding
#' measures of importance see \link[ranger]{ranger}.
#' @param initialRandomEffects Numeric value or vector of initial estimate of random effects.
#' Defaults to 0.
#' @param ErrorTolerance Numeric value to monitor the MERF algorithm's convergence. Defaults to 1e-04.
#' @param MaxIterations Numeric value specifying the maximal amount of iterations for the
#' MERF algorithm. Defaults to 25.
#' @param B Number of bootstrap replications for MSE estimation procedure proposed by
#' Krennmair et al. (202X). Defaults to 100.
#' @param B_adj Number of bootstrap replications for the adjustment of residual variance. Defaults to 100.
#' @param na.rm Logical. Whether missing values should be removed. Defaults to \code{TRUE}.
#' @param ... Additional parameters are directly passed to the random forest \link[ranger]{ranger}.
#' Most important parameters are for instance \code{mtry} (number of variables to possibly split at
#' in each node), or \code{num.trees} (number of trees). For further details on possible parameters
#' see \link[ranger]{ranger} and the example below.
#' @param aggData Logical input indicating whether aggregated covariate information or unit-level covariate information
#' is used. Defaults to \code{FALSE}, assuming unit-level covariate data.
#' @param popnsize data.frame, comprising information of population size of domains.
#' only needed if \code{aggData = TRUE} and a MSE is requested. Please note that the name
#' of the domain identifier must match the column name of \code{smp_data}.
#' @param OOsample_obs Number of Out-of-sample observations taken from the closest area for potentially unsampled
#' areas. Only needed if \code{aggData = TRUE} with default set to 25.
#' @param ADDsamp_obs Number of Out-of-sample observations taken from the closest area if first iteration for the
#' calculation of calibration weights fails. Only needed if \code{aggData = TRUE} with default set to 0.
#' @param w_min Minimal number of covariates from which informative weights are calculated.
#' Only needed if \code{aggData = TRUE}. Defaults to 3.
#'
#' @return An object of class \code{SAEforest} always includes point estimates for disaggregated mean estimates
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
#' The MERF algorithm iteratively optimizes two separate steps: a) the random forest function, assuming the random
#' effects term to be correct and b) estimates the random effects part, assuming the OOB-predictions from the
#' forest to be correct. Overall convergence of the algorithm is monitored by log-likelihood of a joint model
#' of both components. For further details see Krennmair and Schmid or Hajem et. al. (2014).
#'
#' Generally, the MERF requires covariate micro-data. This function, however also allows for the use of aggregated
#' covariate information, by setting \code{aggData = TRUE}. Aggregated covariate information is adaptively
#' incorporated through calibration-weights based on empirical likelihood for the
#' estimation of area-level means. See methodological details in Krennmair et al. (202X).
#'
#' For the estimation of the MSE, residuals are scaled by a bias-corrected residual variance as
#' proposed Krennmair and Schmid (202X). The bootstrap bias correction follows Mendez and Lohr (2011).
#'
#' Note that the \code{MERFmodel} object is a composition of elements from a random forest of class
#' \code{ranger} and a random effects model of class \code{\link[lme4]{merMod}}. Thus, all generic functions are
#' applicable to corresponding objects. For further details on generic functions see \code{\link[ranger]{ranger}}
#' and \code{\link[lme4]{lmer}} as well as the examples below.
#'
#' @references
#' Krennmair, P. and Schmid, T. (202X). WP 1
#'
#' Mendez, G. and Lohr, S. (2011) Paper
#'
#' @seealso \code{\link{SAEforestObject}}, \code{\link[ranger]{ranger}}, \code{\link[lme4]{lmer}}
#' @examples
#' \dontrun{
#' #Loading data
#' data("eusilcA_popAgg")
#' data("eusilcA_pop")
#' data("eusilcA_smp")
#' data("popnsize")
#'
#' income <- eusilcA_smp$eqIncome
#' X_covar <- eusilcA_smp[,-c(1,16,17,18)]
#'
#'#Example 1:
#'#Calculating point estimates and discussing basic generic functions
#'
#' model1 <- SAEforest_mean(Y = income, X = X_covar, dName = "district",
#'                          smp_data = eusilcA_smp, pop_data = eusilcA_pop)
#'
#'#SAEforest generics:
#' summary(model1)
#'
#'
#'#Example 2:
#'#Calculating point + MSE estimates and passing arguments to the forest
#'#Note that B is unrealistically low to improve example speed
#'
#' model2 <- SAEforest_mean(Y = income, X = X_covar, dName = "district",
#'                          smp_data = eusilcA_smp, pop_data = eusilcA_pop,
#'                          MSE = "nonparametric", B = 5, mtry=5,
#'                          num.trees = 100)
#'
#'#SAEforest generics:
#'summary(model2)
#'summarize_indicators(model2, MSE = TRUE, CV =TRUE)
#'
#'#Example 3:
#'#Calculating point + MSE estimates and passing arguments to the forest
#'#Note that B is unrealistically low to improve example speed
#'
#' model3 <- SAEforest_mean(Y = income, X = X_covar, dName = "district",
#'                             smp_data = eusilcA_smp, pop_data = eusilcA_popAgg,
#'                             MSE = "nonparametric", popnsize = popNsize,
#'                             B = 5, mtry=5, num.trees = 100, aggData = TRUE)
#'
#'#SAEforest generics:
#'summary(model3)
#'summarize_indicators(model3, MSE = TRUE, CV =TRUE)
#'
#'}
#'
#' @export

SAEforest_mean <- function(Y, X, dName, smp_data, pop_data, MSE = "none", aggData =FALSE,
                           popnsize = NULL, OOsample_obs = 25, ADDsamp_obs=0, w_min=3,
                           importance = "impurity", initialRandomEffects = 0,
                           ErrorTolerance = 0.0001, MaxIterations = 25,  B=100,
                           B_adj = 100, na.rm =TRUE,...){

  if(aggData == FALSE){

  # ERROR CHECKS OF INPUTS
  #________________________________________

  if(na.rm == TRUE){
    comp_smp <- complete.cases(smp_data)
    smp_data <- smp_data[comp_smp,]
    Y <- Y[comp_smp]
    X <- X[comp_smp,]
    pop_data <- pop_data[complete.cases(pop_data),]
  }

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


  # Point Estimation
  #________________________________________
  mean_preds <- point_mean(Y = Y, X = X, dName = dName, smp_data = smp_data, pop_data = pop_data,
                           initialRandomEffects = initialRandomEffects, ErrorTolerance = ErrorTolerance,
                           MaxIterations = MaxIterations,importance = importance,...)

  data_specs <- sae_specs(dName = dName,cns = pop_data,smp = smp_data)


  if(MSE == "none"){
    result <- list(
      MERFmodel = c(mean_preds[[2]], call = out_call, data_specs = list(data_specs), data=list(smp_data)),
      Indicators = sortAlpha(mean_preds[[1]],dName =dName),
      MSE_Estimates = NULL,
      AdjustedSD = NULL,
      NrCovar = NULL)

    class(result) <- c("SAEforest_mean", "SAEforest")
    return(result)
  }

  # MSE Estimation
  #________________________________________

  if(MSE != "none"){
    print(paste("Error SD Bootstrap started:"))
    adj_SD <- adjust_ErrorSD(Y=Y, X=X, smp_data = smp_data, mod = mean_preds[[2]], B = B_adj, ...)
    print(paste("Bootstrap with", B,"rounds started"))
  }

  if(MSE == "nonparametric"){
    MSE_estims <- MSE_SAEforest_mean_REB(Y=Y, X = X, dName = dName, smp_data = smp_data,
                                         mod=mean_preds[[2]], ADJsd = adj_SD, pop_data = pop_data, B = B,
                                         initialRandomEffects = initialRandomEffects,
                                         ErrorTolerance = ErrorTolerance, MaxIterations = MaxIterations, ...)

    result <- list(
      MERFmodel = c(mean_preds[[2]], call = out_call, data_specs = list(data_specs), data=list(smp_data)),
      Indicators =  sortAlpha(mean_preds[[1]],dName =dName),
      MSE_Estimates =  sortAlpha(MSE_estims, dName =dName),
      AdjustedSD = adj_SD,
      NrCovar = NULL)

    class(result) <- c("SAEforest_mean", "SAEforest")
    return(result)
  }
}

  if(aggData == TRUE){
    # ERROR CHECKS OF INPUTS
    #________________________________________

    if(na.rm == TRUE){
      comp_smp <- complete.cases(smp_data)
      smp_data <- smp_data[comp_smp,]
      Y <- Y[comp_smp]
      X <- X[comp_smp,]
      pop_data <- pop_data[complete.cases(pop_data),]
    }

    input_checks_meanAGG(Y = Y, X = X, dName = dName, smp_data =smp_data, Xpop_agg =pop_data,
                         initialRandomEffects = initialRandomEffects, ErrorTolerance =ErrorTolerance,
                         MaxIterations =MaxIterations, MSE =MSE, B = B, popnsize =popnsize,
                         OOsample_obs =OOsample_obs, ADDsamp_obs = ADDsamp_obs, w_min =w_min,
                         importance = importance, na.rm = na.rm)

    out_call <- match.call()


    # Make domain variable to character and sort data-sets
    smp_data[[dName]] <- factor(smp_data[[dName]], levels=unique(smp_data[[dName]]))
    pop_data[[dName]] <- factor(pop_data[[dName]], levels = unique(pop_data[[dName]]))
    if(MSE != "none"){
      popnsize[[dName]] <- factor(popnsize[[dName]], levels = unique(pop_data[[dName]]))
    }

    # Order Data according to factors to ease MSE estimation
    order_in <- order(smp_data[[dName]])
    smp_data <- smp_data[order_in,]
    X <- X[order_in,]
    Y <- Y[order_in]

    # Point Estimation
    #________________________________________
    mean_preds <- point_meanAGG(Y = Y, X = X, dName = dName, smp_data = smp_data, Xpop_agg = pop_data,
                                   initialRandomEffects = initialRandomEffects, ErrorTolerance = ErrorTolerance,
                                   MaxIterations = MaxIterations, OOsample_obs = OOsample_obs, ADDsamp_obs = ADDsamp_obs, w_min = w_min,
                                   importance = importance, ...)

    data_specs <- sae_specs(dName = dName,cns = pop_data, smp = smp_data)


    if(MSE == "none"){
      result <- list(
        MERFmodel =  c(mean_preds[[2]], call = out_call, data_specs = list(data_specs), data=list(smp_data)),
        Indicators = sortAlpha(mean_preds[[1]],dName =dName),
        MSE_Estimates = NULL,
        AdjustedSD = NULL,
        NrCovar = mean_preds$wAreaInfo)

      class(result) <- c("SAEforest_meanAGG", "SAEforest")
      return(result)
    }

    # MSE Estimation
    #________________________________________

    if(MSE != "none"){
      print(paste("Error SD Bootstrap started:"))
      adj_SD <- adjust_ErrorSD(Y=Y, X=X, smp_data = smp_data, mod = mean_preds[[2]], B = B_adj, ...)
      print(paste("Bootstrap with", B,"rounds started"))
    }

    if(MSE == "nonparametric"){

      MSE_estims <- MSE_SAEforest_aggOOB_wSet(Y=Y, X = X, mod = mean_preds, smp_data = smp_data,
                                              Xpop_agg = pop_data, dName = dName, ADJsd = adj_SD , B=B,
                                              initialRandomEffects = initialRandomEffects,
                                              ErrorTolerance = ErrorTolerance, MaxIterations = MaxIterations,
                                              popnsize = popnsize, ...)

      result <- list(
        MERFmodel =  c(mean_preds[[2]], call = out_call, data_specs = list(data_specs), data=list(smp_data)),
        Indicators = sortAlpha(mean_preds[[1]],dName =dName),
        MSE_Estimates = sortAlpha(MSE_estims,dName =dName) ,
        AdjustedSD = adj_SD,
        NrCovar = mean_preds$wAreaInfo)

      class(result) <- c("SAEforest_meanAGG", "SAEforest")
      return(result)
    }
  }
}






