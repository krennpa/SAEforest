#' Main function for unit-level MERF
#'
#' This function enables the use of Mixed Effects Random Forests (MERFs) by effectively
#' combining a random forest from \pkg{ranger} with a model capturing random effects from
#' \pkg{lme4}. The MERF algorithm is an algorithmic procedure reminiscent of an EM-algorithm
#' (see Details). The function is the base-function for the wrapping function (\code{\link{SAEforest_model}}
#' and should not be directly used by the ordinary user. Recommended exceptions are applications exceeding
#' the scope of existing wrapper functions or further research. The function \code{MERFranger}
#' allows to model complex patterns of structural relations (see Examples). The function returns
#' an object of class \code{MERFranger}, which can be used to produce unit-level predictions. In contrast to
#' the wrapping functions, this function does not directly provide SAE estimates on domain-specific indicators.
#'
#' @param Y Continuous input value of target variable.
#' @param X Matrix of predictive covariates.
#' @param random Specification of random effects terms following the syntax of \link[lme4]{lmer}.
#' Random effect terms are specified by vertical bars \code{(|)} separating expressions for design matrices
#' from grouping factors. For further details see \link[lme4]{lmer} and the example below.
#' @param data data.frame of sample data including the specified elements of \code{Y} and
#' \code{X}.
#' @param initialRandomEffects Numeric value or vector of initial estimate of random effects.
#' Defaults to 0.
#' @param ErrorTolerance Numeric value to monitor the MERF algorithm's convergence. Defaults to 1e-04.
#' @param MaxIterations Numeric value specifying the maximal amount of iterations for the
#' MERF algorithm. Defaults to 25.
#' @param importance Variable importance mode processed by the
#' random forest from the \pkg{ranger}. Must be 'none', 'impurity', 'impurity_corrected',
#' 'permutation'. For further details see \link[ranger]{ranger}.
#' @param na.rm Logical. Whether missing values should be removed. Defaults to \code{TRUE}.
#' @param ... Additional parameters are directly passed to the random forest \link[ranger]{ranger}.
#' Most important parameters are for instance \code{mtry} (number of variables to possibly split at
#' in each node), or \code{num.trees} (number of trees). For further details on possible parameters
#' see \link[ranger]{ranger} and the example below.
#'
#' @return An object of class MERFranger includes the following elements:
#'
#' \item{\code{Forest}}{A random forest of class \link[ranger]{ranger} modelling fixed effects
#' of the model.}
#' \item{\code{EffectModel}}{A model of random effects of class \code{\link[lme4]{merMod}} capturing
#' structural components of MERFs and modeling random components.}
#' \item{\code{RandomEffects}}{List element containing the values of random intercepts from \code{EffectModel}.}
#' \item{\code{RanEffSD}}{Numeric value of the standard deviation of random intercepts.}
#' \item{\code{ErrorSD}}{Numeric value of standard deviation of unit-level errors.}
#' \item{\code{VarianceCovariance}}{VarCorr matrix from \code{EffectModel}.}
#' \item{\code{LogLik}}{Vector with numerical entries showing the loglikelihood of the MERF algorithm.}
#' \item{\code{IterationsUsed}}{Numeric number of iterations used until convergence of the MERF algorithm.}
#' \item{\code{OOBresiduals}}{Vector of OOB-residuals.}
#' \item{\code{Random}}{Character specifying the random intercept in the random effects model.}
#' \item{\code{ErrorTolerance}}{Numerical value to monitor the MERF algorithm's convergence.}
#' \item{\code{initialRandomEffects}}{Numeric value or vector of initial specification of random effects.}
#' \item{\code{MaxIterations}}{Numeric value specifying the maximal amount of iterations for the
#' MERF algorithm.}
#'
#' @details
#' There exists a generic function for \code{predict} for objects obtained by \code{MERFranger}.
#'
#' The MERF algorithm iteratively optimizes two separate steps: a) the random forest
#' function, assuming the random effects term to be correct and b) estimates the random
#' effects part, assuming the OOB-predictions from the forest to be correct. Overall convergence
#' of the algorithm is monitored by the log-likelihood of a joint model of both components. For
#' further details see Krennmair & Schmid (2022) or Hajjem et al. (2014).
#'
#' Note that the \code{MERFranger} object is a composition of elements from a random forest of class
#' \code{ranger} and a random effects model of class \code{\link[lme4]{merMod}}. Thus, all generic functions are
#' applicable to corresponding objects. For further details on generic functions see \code{\link[ranger]{ranger}}
#' and \code{\link[lme4]{lmer}} as well as the examples below.
#'
#' @references
#' Hajjem, A., Bellavance, F., & Larocque, D. (2014). Mixed-Effects Random Forest for Clustered
#' Data. Journal of Statistical Computation and Simulation, 84 (6), 1313–1328.
#'
#' Krennmair, P., & Schmid, T. (2022). Flexible Domain Prediction Using Mixed Effects
#' Random Forests. Journal of Royal Statistical Society: Series C (Applied Statistics) (forthcoming).
#'
#' @seealso \code{\link{SAEforest}}, \code{\link[ranger]{ranger}}, \code{\link[lme4]{lmer}},
#' \code{\link{SAEforest_model}}
#'
#' @export
#' @import ranger
#' @import lme4
#'
#' @examples
#' \donttest{
#' # Load Data
#' data("eusilcA_pop")
#' data("eusilcA_smp")
#'
#' income <- eusilcA_smp$eqIncome
#' X_covar <- eusilcA_smp[, -c(1, 16, 17, 18)]
#'
#' # Example 1:
#' # Calculating general model used in wrapper functions
#'
#' model1 <- MERFranger(Y = income, X = X_covar, random = "(1|district)",
#'                      data = eusilcA_smp, num.trees=50)
#'
#' # get individual predictions:
#'
#' ind_pred <- predict(model1, eusilcA_pop)
#'}
MERFranger <- function(Y,
                       X,
                       random,
                       data,
                       importance = "none",
                       initialRandomEffects = 0,
                       ErrorTolerance = 0.0001,
                       MaxIterations = 25,
                       na.rm = TRUE,
                       ...) {

  if (na.rm == TRUE) {
    comp_smp <- complete.cases(data)
    data <- data[comp_smp, ]
    Y <- Y[comp_smp]
    X <- X[comp_smp, ]
  }

  input_checks_MERF(
    Y = Y, X = X, data = data, initialRandomEffects = initialRandomEffects,
    ErrorTolerance = ErrorTolerance, MaxIterations = MaxIterations,
    importance = importance, na.rm = na.rm
  )

  Target <- Y
  ContinueCondition <- TRUE
  iterations <- 0

  AdjustedTarget <- Target - initialRandomEffects
  oldLogLik <- 0
  while (ContinueCondition) {
    iterations <- iterations + 1
    rf <- ranger::ranger(x = X, y = AdjustedTarget, importance = importance, ...)
    forest_preds <- rf$predictions
    f0 <- as.formula(paste0("Target ~ -1+", random))
    lmefit <- lme4::lmer(f0, data = data, REML = FALSE, offset = forest_preds)

    newLogLik <- as.numeric(stats::logLik(lmefit))

    ContinueCondition <- (abs(abs(newLogLik - oldLogLik[iterations]) / oldLogLik[iterations]) > ErrorTolerance &
      iterations < MaxIterations)
    oldLogLik <- c(oldLogLik, newLogLik)
    AllEffects <- predict(lmefit)
    AdjustedTarget <- Target - (AllEffects - forest_preds)
  }

  data$forest_preds <- NULL
  residuals <- Target - predict(lmefit)

  result <- list(
    Forest = rf,
    EffectModel = lmefit,
    RandomEffects = lme4::ranef(lmefit),
    RanEffSD = as.data.frame(lme4::VarCorr(lmefit))$sdcor[1],
    ErrorSD = stats::sigma(lmefit),
    VarianceCovariance = lme4::VarCorr(lmefit),
    LogLik = oldLogLik,
    IterationsUsed = iterations,
    OOBresiduals = residuals,
    Random = random,
    ErrorTolerance = ErrorTolerance,
    initialRandomEffects = initialRandomEffects,
    MaxIterations = MaxIterations
  )

  class(result) <- "MERFranger"

  return(result)
}

# Generic predict function for MERFranger ------------------------------------------------
#' @export
predict.MERFranger <- function(object, ...) {
  retval <- predict(object$Forest, ...)$predictions +
    predict(object$EffectModel, allow.new.levels = TRUE, ...)

  return(retval)
}
