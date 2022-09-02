#' Presents point, MSE and CV estimates
#'
#' Function \code{summarize_indicators} reports point and
#' mean squared error (MSE) estimates as well as calculated coefficients of variation
#' (CV) from a fitted \code{SAEforest} object.
#'
#' @param object Object for which point and/or MSE estimates and/or
#' calculated CV's are requested. The object must be of class \code{SAEforest}.
#' @param indicator Optional character vector specifying indicators to be mapped: (i)
#' all calculated indicators ("all"); (ii) each default indicators name: "Mean",
#' "Quant10", "Quant25", "Median", "Quant75", "Quant90", "Gini", "Hcr", "Pgap", "Qsr"
#' or the function name/s of "custom_indicator/s"; (iii) a vector of names of indicators.
#' If the \code{object} is estimated by \code{\link{SAEforest_model}} indicator arguments
#' are ignored and only the "Mean" is returned.
#' @param MSE Logical. If \code{TRUE}, MSE estimates for selected indicators
#' per domain are added to the data frame of point estimates. Defaults to \code{FALSE}.
#' @param CV Logical. If \code{TRUE}, coefficients of variation for selected
#' indicators per domain are added to the data frame of point estimates. Defaults to \code{FALSE}.
#'
#' @return The return of \code{summarize_indicators} is an object of class \code{summarize_indicators.SAEforest}
#' including domain-specific point and/or MSE estimates and/or calculated CV's from a \code{SAEforest} object
#' The returned object contains the data.frame \code{ind} and a character including the names of requested indicator(s).
#'
#' @details Objects of class \code{summarize_indicators.SAEforest} have methods for following generic
#' functions: \code{head} and \code{tail} (for default documentation, see
#' \code{\link[utils]{head}}),  \code{as.matrix} (for default documentation, see
#' \code{\link[base]{matrix}}), \code{as.data.frame} (for default documentation,
#' see \code{\link[base]{as.data.frame}}), \code{subset} (for default
#' documentation, see \code{\link[base]{subset}}).
#'
#' @seealso \code{\link{SAEforestObject}}, \code{\link{SAEforest_model}}
#'
#' @examples
#' \donttest{
#' # Loading data
#' data("eusilcA_pop")
#' data("eusilcA_smp")
#'
#' income <- eusilcA_smp$eqIncome
#' X_covar <- eusilcA_smp[, -c(1, 16, 17, 18)]
#'
#' # Calculating point + MSE estimates and passing arguments to the forest.
#' # Additionally, two additional indicators and functions as threshold are added.
#' # Note that B and num.trees are low to speed up estimation time and must be changed for
#' # practical applications.
#'
#' model1 <- SAEforest_model(Y = income, X = X_covar, dName = "district",
#'                           smp_data = eusilcA_smp, pop_data = eusilcA_pop,
#'                           meanOnly = FALSE, MSE = "nonparametric", B = 5, mtry = 5,
#'                           num.trees = 50, smearing = FALSE)
#'
#' # Extract indicator and show generics:
#' Gini1 <- summarize_indicators(model1, MSE = TRUE, CV = TRUE, indicator = "Gini")
#'
#' head(Gini1)
#' tail(Gini1)
#' as.data.frame(Gini1)
#' as.matrix(Gini1)
#' subset(Gini1, district == "Wien")
#' }
#' @export


summarize_indicators <- function(object,
                                 indicator = "all",
                                 MSE = FALSE,
                                 CV = FALSE) {

  summarize_indicators_check(object = object, indicator = indicator, MSE = MSE, CV = CV)

  if ( (sum(indicator == "all")>=1) || sum((indicator == "All")>=1)) {
    indicator <- colnames(object$Indicators)[-1]
  }

  if (inherits(object, "SAEforest_mean") || inherits(object, "SAEforest_meanAGG")) {
    indicator <- "Mean"
  }

  out_var <- data.frame(
    object$Indicators[colnames(object$Indicators)[1]],
    object$Indicators[indicator]
  )

  selected <- colnames(out_var)[-1]

  if (MSE == TRUE || CV == TRUE) {
    MSE_estims <- object$MSE_Estimates[indicator]
    cv_estims <- sqrt(object$MSE_Estimates[indicator]) / object$Indicators[indicator]
    colnames(cv_estims) <- indicator
    colnames(MSE_estims) <- paste0(colnames(MSE_estims), "_MSE")
    colnames(cv_estims) <- paste0(colnames(cv_estims), "_CV")
    combined <- data.frame(out_var, MSE_estims, cv_estims)
    endings <- c("", "_MSE", "_CV")[c(TRUE, MSE, CV)]

    combined <- combined[, c(paste0(colnames(combined)[1]), paste0(
      rep(selected, each = length(endings)),
      endings
    ))]
  } else {
    combined <- out_var
  }
  estimators_SAEforest <- list(ind = combined, ind_name = indicator)

  class(estimators_SAEforest) <- "summarize_indicators.SAEforest"
  return(estimators_SAEforest)
}


# Generic Functions for summarize_indicators ----------------------------------------------

# Tail/head functions
# Prints summarize_indicators.SAEforest objects
#' @export

print.summarize_indicators.SAEforest <- function(x, ...) {
  cat(paste0("Indicator/s: ", x$ind_name, "\n"))
  print(x$ind)
}


#' @importFrom utils head
#' @export
# CV estimators

head.summarize_indicators.SAEforest <- function(x, n = 6L, addrownums = NULL, ...) {
  head(x$ind, n = n, addrownums = addrownums, ...)
}

#' @importFrom utils tail
#' @export

tail.summarize_indicators.SAEforest <- function(x, n = 6L, keepnums = TRUE, addrownums = NULL, ...) {
  tail(x$ind, n = n, keepnums = keepnums, ...)
}

# Transforms summarize_indicators.SAEforest objects into a matrix object
#' @export

as.matrix.summarize_indicators.SAEforest <- function(x, ...) {
  as.matrix(x$ind[, -1])
}

# Transforms summarize_indicators.SAEforest objects into a dataframe object
#' @export

as.data.frame.summarize_indicators.SAEforest <- function(x, ...) {
  as.data.frame(x$ind, ...)
}

# Subsets an summarize_indicators.SAEforest object
#' @export

subset.summarize_indicators.SAEforest <- function(x, ...) {
  x <- as.data.frame(x)
  subset(x = x, ...)
}
