# Definition of some useful generics from lme4

# Extract fixed effects of a SAEforest object ---------------------------------------------
# Output corresponds to an object of class ranger
#' @export
#' @method fixef SAEforest

fixef.SAEforest <- function(object, add.dropped = FALSE, ...) {
  class_error(object)
  object$MERFmodel$Forest
}


# Extract survey data of a SAEforest object ----------------------------------------------
#' @export
#' @method getData SAEforest

getData.SAEforest <- function(object) {
  class_error(object)
  object$MERFmodel$data
}


# Extract ranef of a SAEforest object ----------------------------------------------------
#' @export
#' @method ranef SAEforest

ranef.SAEforest <- function(object, ...) {
  class_error(object)
  ranef(object$MERFmodel$EffectModel, ...)
}


# Extract residuals of a SAEforest object ------------------------------------------------
#' @export
#' @method residuals SAEforest

residuals.SAEforest <- function(object, ...) {
  class_error(object)
  object$MERFmodel$OOBresiduals
}


# Extract sigma of a SAEforest object ----------------------------------------------------
#' @importFrom stats sigma
#' @export
#' @method sigma SAEforest

sigma.SAEforest <- function(object, ...) {
  class_error(object)
  sigma(object$MERFmodel$EffectModel, ...)
}


# Extract VarCorr of a SAEforest object --------------------------------------------------
#' @export
#' @method VarCorr SAEforest

VarCorr.SAEforest <- function(x, sigma = 1, ...) {
  class_error(x)
  VarCorr(x$MERFmodel$EffectModel, ...)
}
