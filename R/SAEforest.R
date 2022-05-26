#' SAEforest - Estimating disaggregated indicators using Mixed Effects Random Forests
#'
#' The package \pkg{SAEforest} promotes the use of Mixed Effects Random Forests (MERFs) for
#' applications of Small Area Estimation (SAE). The package effectively combines functions for
#' the estimation of regionally disaggregated linear and non-linear economic and inequality
#' indicators using survey sample data. Estimated models increase the precision of direct
#' estimates from survey data, combining unit-level and aggregated population level covariate
#' information from census or register data. Apart from point estimates, MSE estimates for requested
#' indicators can be easily obtained. The package provides procedures to facilitate the analysis
#' of model performance of MERFs and visualizes predictive relations from covariates and variable
#' importance. Additionally, users can summarize and map indicators and corresponding measures of
#' uncertainty.
#'
#' @details
#' This package includes a main function \code{\link{MERFranger}} that is wrapped into two functions
#' \code{\link{SAEforest_mean}} and \code{\link{SAEforest_nonLin}} for an improved SAE workflow.
#' Each function produces an object inheriting requested results of regionally disaggregated point
#' and uncertainty estimates. Additionally, statistical information on model fit and variable
#' importance is accessible through generic functions such as a summary (\code{\link{summary.SAEforest}})
#' or a class specific plot function (\code{\link{plot.SAEforest}}). For a full documentation of
#' objects of class \code{SAEforest} see \code{\link{SAEforestObject}}. An overview of all currently
#' provided functions within this package can be be seen with \code{help(package="SAEforest")}.
#'
#' @references
#' Krennmair, P. and Schmid, T. (202X). WP 1
#' Krennmair, P. Würz, N. and Schmid, T. (202X). WP 2
#' Krennmair, P., Schmid, T. and Tzavidis N. (202X). WP 3
#'
#' @docType package
#' @name SAEforest
NULL
