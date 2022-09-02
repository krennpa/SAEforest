#' Prints a 'SAEforest' object
#'
#' Basic information of an \code{\link{SAEforestObject}} is printed.
#' @param x Object of class \code{SAEforest}, representing point and MSE
#' estimates obtained by function \code{\link{SAEforest_model}}.
#' @seealso \code{\link{SAEforestObject}}
#' @param ... Optional additional inputs that are ignored for this method.
#'
#' @return Prints basic information on survey data characteristics.
#'
#' @export
print.SAEforest <- function(x, ...) {
  class_error(x)
  cat("________________________________________________________________\n")
  cat("Mixed Effects Random Forest for Small Area Estimation\n")
  cat("________________________________________________________________\n")
  cat("\n")
  cat("Information on Domains\n")
  total_dom <- x$MERFmodel$data_specs$D_total
  in_dom <- x$MERFmodel$data_specs$D_in
  oos_dom <- x$MERFmodel$data_specs$D_out

  dom_info <- data.frame(in_dom, oos_dom, total_dom)
  rownames(dom_info) <- c("")
  colnames(dom_info) <- c("In-sample", "Out-of-sample", "Total")

  smp_size <- x$MERFmodel$data_specs$N_surv
  pop_size <- x$MERFmodel$data_specs$N_pop

  print(dom_info)
  cat("\n")
  cat("Units in sample:", smp_size, "\n")
  if (!inherits(x, "SAEforest_meanAGG")) {
    cat("Units in population:", pop_size, "\n")
  }
  cat("\n")
  cat("Model Information:\n")
  cat("An SAEforest Object contains results for point- and uncertainty estiamtes as well as seperate\n")
  cat("model components of the random forest part as well as the mixed effects part\n")
  cat("\n")
  cat("Mehtods of lme4 are applicable to the random effects components stored in 'object$MERFmodel$EffectModel'.\n")
  cat("Mehtods for random forests from ranger are applicable to the fixed effects components stored in\n")
  cat("'object$MERFmodel$Forest'.\n")
}
