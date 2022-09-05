#' Summarizes an 'SAEforest' object
#'
#' Shows additional information about the data, the SAE model and its components.
#' Information is extracted from a \code{SAEforest} object. The returned object
#' is suitable for printing with \code{print}.

#' @param object An object of class \code{SAEforest} representing point
#' and MSE estimates. Objects differ depending on the estimation method.
#' @param ... Optional additional inputs that are ignored for this method.
#'
#' @return An object of class \code{summary.SAEforest} including information about the sample
#' and population data, the model fit and random forest specific metrics.
#' @seealso \code{\link{SAEforestObject}}
#' @examples
#' \donttest{
#' # Loading data
#' data("eusilcA_pop")
#' data("eusilcA_smp")
#'
#' income <- eusilcA_smp$eqIncome
#' X_covar <- eusilcA_smp[, -c(1, 16, 17, 18)]
#'
#' # Example 1:
#' # Calculating point estimates and discussing basic generic functions
#'
#' model1 <- SAEforest_model(Y = income, X = X_covar, dName = "district",
#'                          smp_data = eusilcA_smp, pop_data = eusilcA_pop,
#'                          num.trees=50, mtry=3)
#'
#' # SAEforest generics:
#' summary(model1)
#' }
#' @export

summary.SAEforest <- function(object, ...) {
  class_error(object = object)

  call_SAE <- object$MERFmodel$call

  total_dom <- object$MERFmodel$data_specs$D_total
  in_dom <- object$MERFmodel$data_specs$D_in
  oos_dom <- object$MERFmodel$data_specs$D_out

  dom_info <- data.frame(in_dom, oos_dom, total_dom)
  rownames(dom_info) <- c("")
  colnames(dom_info) <- c("In-sample", "Out-of-sample", "Total")

  smp_size <- object$MERFmodel$data_specs$N_surv
  pop_size <- object$MERFmodel$data_specs$N_pop

  if (inherits(object, "SAEforest_meanAGG")) {
    pop_size <- NULL
  }

  smp_size_dom <- summary(as.data.frame(object$MERFmodel$data_specs$ni_smp)[, "Freq"])
  pop_size_dom <- summary(as.data.frame(object$MERFmodel$data_specs$ni_pop)[, "Freq"])

  sizedom_smp_pop <- rbind(
    Sample_domains = smp_size_dom,
    Population_domains = pop_size_dom
  )

  if (inherits(object, "SAEforest_meanAGG")) {
    sizedom_smp_pop <- rbind(Sample_domains = smp_size_dom)
  }

  # information on forest:
  forest_info <- data.frame(c(
    object$MERFmodel$Forest$treetype, object$MERFmodel$Forest$num.trees, object$MERFmodel$Forest$num.independent.variables,
    object$MERFmodel$Forest$mtry, object$MERFmodel$Forest$min.node.size, object$MERFmodel$Forest$importance.mode,
    object$MERFmodel$Forest$splitrule, round(object$MERFmodel$Forest$r.squared, digits = 5)
  ))

  colnames(forest_info) <- NULL
  rownames(forest_info) <- c(
    "Type:", "Number of trees:", "Number of independent variables:", "Mtry:",
    "Minimal node size:", "Variable importance mode:", "Splitrule:", "Rsquared (OOB):"
  )

  # information on LM:
  sum_lm <- summary(object$MERFmodel$EffectModel)

  # ICC
  icc <- object$MERFmodel$RanEffSD^2 / (object$MERFmodel$RanEffSD^2 + object$MERFmodel$ErrorSD^2)

  # information on convergence:
  LogLik <- as.data.frame(t(object$MERFmodel$LogLik))
  rownames(LogLik) <- c("")
  colnames(LogLik) <- NULL
  iter <- object$MERFmodel$IterationsUsed
  maxIter <- object$MERFmodel$MaxIterations
  Tol <- object$MERFmodel$ErrorTolerance

  sum_SAEforest <- list(
    call_SAE = call_SAE,
    dom_info = dom_info,
    smp_size = smp_size,
    pop_size = pop_size,
    sizedom_smp_pop = sizedom_smp_pop,
    forest_info = forest_info,
    sum_lm = sum_lm,
    icc = icc,
    LogLik = LogLik,
    iter = iter,
    maxIter = maxIter,
    Tol = Tol
  )

  class(sum_SAEforest) <- c("summary.SAEforest", "SAEforest")
  sum_SAEforest
}

# Generic print function for summary.SAEforest --------------------------------------------

#' @export
print.summary.SAEforest <- function(x, ...) {
  class_error(object = x)
  cat("________________________________________________________________\n")
  cat("Mixed Effects Random Forest for Small Area Estimation\n")
  cat("________________________________________________________________\n")
  cat("Call:\n")
  print(x$call_SAE)
  cat("\n")
  cat("Domains\n")
  cat("________________________________________________________________")
  cat("\n")
  print(x$dom_info)
  cat("\n")
  cat("Totals:\n")
  cat("Units in sample:", x$smp_size, "\n")
  if (!is.null(x$pop_size)) {
    cat("Units in population:", x$pop_size, "\n")
  }
  cat("\n")
  print(x$sizedom_smp_pop)
  cat("\n")
  cat("Random forest component: \n")
  cat("________________________________________________________________\n")
  print(x$forest_info)
  cat("\n")
  cat("Structural component of random effects:\n")
  cat("________________________________________________________________\n")
  print(x$sum_lm)
  cat("\n")
  cat("ICC: ", x$icc, "\n")
  cat("\n")
  cat("Convergence of MERF algorithm: \n")
  cat("________________________________________________________________\n")
  cat("Convergence achieved after", x$iter, "iterations.\nA maximum of", x$maxIter, "iterations used and tolerance set to:", x$Tol, "\n")
  cat("\n")
  cat("Monitored Log-Likelihood:")
  print(x$LogLik)
}
