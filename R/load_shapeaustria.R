#' Loading the shape file for Austrian districts
#'
#' The function simplifies the loading of the shape file for Austrian districts.
#' It is originally used for examples in package \pkg{emdi}.
#'
#' @return A shape file of class \code{SpatialPolygonsDataFrame}.
#' @details The shape file contains the borders of 94 Austrian districts.
#' The main purpose of this function is the visualization of estimation results with the
#' plotting function \code{\link{map_indicators}}.
#' @seealso Information on the class of \code{\link[sp]{SpatialPolygonsDataFrame}} from the package \pkg{sp}.
#'
#' @export

load_shapeaustria <- function(){
  load(system.file("shapes/shape_austria_dis.rda", package = "SAEforest"),
       envir = .GlobalEnv)
}
