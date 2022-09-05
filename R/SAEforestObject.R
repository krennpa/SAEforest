#' Fitted 'SAEforest' object
#'
#' An object of class \code{SAEforest} always includes point estimates of regionally disaggregated economic
#' and inequality indicators and a \code{MERFmodel} element including information on the model fit for fixed
#' effects as well as random effects. Optionally an \code{SAEforestObject} includes corresponding MSE estimates.
#' In the case of mean estimates and aggregated covariate information, the \code{SAEforestObject} additionally
#' includes an element, capturing the number of variables used in the weighting process from aggregated
#' covariate information. For an object of class \code{SAEforestObject}, the following generic functions are applicable:
#' \code{\link{print}}, \code{\link{plot}}, \code{\link{summary}} and \code{\link{summarize_indicators}}.
#' Additionally selected generic functions of \pkg{lme4} (\code{fixef, getData, ranef, residuals, sigma, VarCorr}) are
#' directly applicable to an object of class \code{SAEforest}.
#'
#' @return
#' Four components are always included in an SAEforest object. \code{MSE_estimates} and \code{AdjustedSD} are
#' \code{NULL} except MSE results are requested. An element of \code{NrCovar} only exists for SAEforest objects
#' produced by \code{\link{SAEforest_model}} with option \code{aggData = TRUE}.
#'
#' \item{\code{MERFmodel}}{The included \code{MERFmodel} object comprises information on the model fit, details
#' on the performed MERF algorithm as well as details on variance components. See below for an exact description
#' of components.}
#' \item{\code{Indicators}}{A data frame where the first column is the area-level identifier and additional columns
#' are the indicators of interest. Note that objects from \code{\link{SAEforest_model}}
#' only report the "Mean".}
#' \item{\code{MSE_estimates}}{Only if MSE results requested. A data frame where the first column is the area-level
#' identifier and additional columns are the MSE estimates for indicators of interest. Note that objects from
#' \code{\link{SAEforest_model}} only report MSE values for the "Mean".}
#' \item{\code{NrCovar}}{Only if means under aggregated covariate information are estimated, i.e.
#' \code{\link{SAEforest_model}} with option \code{aggData = TRUE}. A list containing variable names of
#' covariates used for the calculation of needed calibration weights for point estimates. See Krennmair et al. (2022a) for
#' methodological details an explanations.}
#'
#' Details on object of \code{MERFmodel}:
#'
#' \item{\code{Forest}}{A random forest of class \link[ranger]{ranger} modelling fixed effects
#' of the model.}
#' \item{\code{EffectModel}}{A model of random effects of class \code{\link[lme4]{merMod}} capturing
#' structural components of MERFs and modeling random components.}
#' \item{\code{RandomEffects}}{List element containing the values of random intercepts from \code{EffectModel}.}
#' \item{\code{RanEffSD}}{Numeric value of standard deviation of random intercepts.}
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
#' \item{\code{call}}{The summarized function call producing the object.}
#' \item{\code{data_specs}}{Data characteristics such as domain-specific sample sizes or number of
#' out-of-sample areas.}
#' \item{\code{data}}{Processed survey sample data.}
#'
#' @details
#' Note that the \code{MERFmodel} object is a composition of elements from a random forest of class
#' \code{ranger} and a random effects model of class \code{\link[lme4]{merMod}}. Thus, all generic functions are
#' applicable to corresponding objects. For further details on generic functions see \code{\link[ranger]{ranger}}
#' and \code{\link[lme4]{lmer}} as well as the examples below.
#'
#' @references
#' Krennmair, P., & Schmid, T. (2022). Flexible Domain Prediction Using Mixed Effects
#' Random Forests. Journal of Royal Statistical Society: Series C (Applied Statistics) (forthcoming).
#'
#' Krennmair, P., & Würz, N. & Schmid, T. (2022a). Analysing Opportunity Cost of Care Work using
#' Mixed Effects Random Forests under Aggregated Census Data.
#'
#' Krennmair, P., & Schmid, T & Tzavidis, Nikos. (2022b). The Estimation of Poverty Indicators Using
#' Mixed Effects Random Forests. Working Paper.
#'
#' @seealso \code{\link{SAEforest_model}}, \code{ \link[ranger]{ranger}},
#' \code{ \link[lme4]{lmer}}
#'
#' @examples
#' \donttest{
#' # Loading data
#' data("eusilcA_pop")
#' data("eusilcA_smp")
#'
#' income <- eusilcA_smp$eqIncome
#' X_covar <- eusilcA_smp[,-c(1,16,17,18)]
#'
#' # Example 1:
#' # Calculating point estimates and discussing basic generic functions
#'
#' model1 <- SAEforest_model(Y = income, X = X_covar, dName = "district",
#'                          smp_data = eusilcA_smp, pop_data = eusilcA_pop,
#'                          num.trees=50, mtry = 3)
#'
#' #SAEforest generics:
#'
#' summary(model1)
#' summarize_indicators(model1)
#' residuals(model1)
#' sigma(model1)
#' }
#' @name SAEforestObject
NULL
