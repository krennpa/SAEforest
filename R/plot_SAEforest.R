#' Plot function for a 'SAEforest' object
#'
#' Plots model-specific characteristics of the fixed effects random forest component of
#' the MERF from a \code{\link{SAEforestObject}}. A variable importance plot is produced to visualize
#' the importance of individual covariates for the predictive performance of the model.
#' For the variable importance plot, arguments are passed internally to the function
#' \code{\link[vip]{vip}}. If requested, the plot function additionally provides a partial
#' dependence plot (pdp) to visualize the impact of a given number of influential covariates
#' on the target variable. The pdp plot is produced using \code{\link[pdp]{partial}} from
#' the package \pkg{pdp}. The plot-engine for both plots is \pkg{ggplot2}.
#'
#' @param x An object of class \code{SAEforest} including a random forest model of class \code{\link[ranger]{ranger}}.
#' @param num_features Number of features for which a partial dependence plot is required.
#' @param col Parameter specifying the color of selected plots. The argument must be specified
#' such that it can be processed by \code{\link[ggplot2]{aes}}. Defaults to a character name of the
#' color "darkgreen".
#' @param fill Parameter specifying the fill of selected plots. The argument must be specified
#' such that it can be processed by \code{\link[ggplot2]{aes}}. Defaults to a character name of the
#' color "darkgreen".
#' @param alpha Parameter specifying the transparency of \code{fill} for \code{\link[vip]{vip}} plots.
#' The argument must be a number in \code{[0,1]}.
#' @param include_type Logical. If set to \code{TRUE}, the type of importance specified in the fitting process
#' of the model is included in the \code{\link[vip]{vip}} plot. Defaults to \code{TRUE}.
#' @param horizontal Logical. If set to \code{TRUE}, the importance scores appear on the x-axis. If parameter is
#' set to \code{FALSE}, the importance scores are plot on the y-axis. Defaults to \code{TRUE}.
#' @param gg_theme Specify a predefined theme from \pkg{ggplot2}. Defaults to \code{theme_minimal}.
#' @param lsize Parameter specifying the line size of pdp plots. The argument must be specified
#' such that it can be processed by \code{\link[ggplot2]{aes}}. Defaults to 1.5.
#' @param lty Parameter specifying the line size of pdp plots. The argument must be specified
#' such that it can be processed by \code{\link[ggplot2]{aes}}. Defaults to "solid".
#' @param grid_row Parameter specifying the amount of rows for the joint pdp plot. Defaults to 2.
#' @param out_list Logical. If set to \code{TRUE}, a list of individual plots produced by \pkg{ggplot2}
#' is returned for further individual customization and processing. Defaults to \code{FALSE}.
#' @param pdp_plot Logical. If set to \code{TRUE}, partial dependence plots produced by \code{\link[pdp]{partial}}
#' from the package \pkg{pdp} are included. Defaults to \code{TRUE}.
#' @param ... Optional additional inputs that are ignored for this method.
#'
#' @return Plots of variable importance and/or partial dependence of covariates ranked by corresponding
#' importance. Additionally, a list of individual plots can be returned facilitating individual
#' customization and exporting. See the following examples for details.
#'
#' @details For the production of importance plots, be sure to specify the parameter of
#' \code{importance != 'none'} before producing estimates with function \code{\link{SAEforest_model}}.
#'
#' For pdp plots, note that covariates of type factor or character cannot be used for partial dependence
#' plots. Dummy-variables can be used, however, their pdp plots are always lines connecting two effect
#' points for 0 and 1. Most informative pdp plots can be produced for continuous predictors.
#'
#' @seealso \code{\link{SAEforestObject}}
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
#' # Example 1:
#' # Calculating point estimates and discussing basic generic functions
#'
#' model1 <- SAEforest_model(Y = income, X = X_covar, dName = "district",
#'                           smp_data = eusilcA_smp, pop_data = eusilcA_pop,
#'                           num.trees = 50)
#' plot(model1)
#' }
#'
#' @importFrom ggplot2 theme_minimal ggplot aes_string ggtitle geom_line
#' @export
plot.SAEforest <- function(x,
                           num_features = 6,
                           col = "darkgreen",
                           fill = "darkgreen",
                           alpha = 0.8,
                           include_type = TRUE,
                           horizontal = TRUE,
                           gg_theme = theme_minimal(),
                           lsize = 1.5,
                           lty = "solid",
                           grid_row = 2,
                           out_list = FALSE,
                           pdp_plot = TRUE,
                           ...) {
  class_error(x)

  input_checks_plot(
    num_features = num_features, alpha = alpha, include_type = include_type, horizontal = horizontal,
    lsize = lsize, grid_row = grid_row, out_list = out_list, pdp_plot = pdp_plot, gg_theme = gg_theme
  )

  # produce a vip plot
  vip_plot <- vip::vip(x$MERFmodel$Forest,
                       aes = list(col = col, fill = fill, alpha = alpha),
                       include_type = include_type, horizontal = horizontal, num_features = num_features
  ) + ggtitle("Variable Importance") + gg_theme

  eval(print(vip_plot))
  cat("Press [enter] to continue")
  line <- readline()

  # produce a pdp plot and check if variables are factors or characters
  pdp_curves <- NULL

  if (pdp_plot == TRUE) {
    set_fact <- names(x$MERFmodel$data)[sapply(x$MERFmodel$data, is.factor)]
    set_char <- names(x$MERFmodel$data)[sapply(x$MERFmodel$data, is.character)]

    set_rm <- levels(factor(c(set_fact, set_char)))

    if (length(set_rm) != 0) {
      warning(paste0("The data contained ", length(set_rm), " character or factor variables unsuitable for pdp plots(", paste(set_rm, collapse = ", "), ")."))
    }

    forest_imp <- as.data.frame(vip::vi(x$MERFmodel$Forest))
    forest_imp <- forest_imp[order(forest_imp$Importance, decreasing = TRUE), ]
    forest_imp <- forest_imp[!forest_imp[, "Variable"] %in% set_rm, ]
    forest_imp <- na.omit(forest_imp[1:num_features, ])
    y_name <- as.character(x$MERFmodel$call$Y)
    y_name <- paste0(gsub("[[:punct:]]", "", y_name),collapse="")

    rm_extraLabels <- haven::zap_labels(x$MERFmodel$data)

    pdp_curves <- lapply(forest_imp[, "Variable"], FUN = function(feature) {
      pd <- pdp::partial(x$MERFmodel$Forest, pred.var = feature, train = rm_extraLabels, plot = FALSE)
      colnames(pd)[2] <- y_name
      ggplot(data = pd, aes_string(y = y_name, x = feature)) +
        geom_line(linetype = lty, color = col, size = lsize) +
        ggtitle(paste("Partial Dependence of", feature)) +
        gg_theme
    })

    vip::grid.arrange(grobs = pdp_curves, nrow = grid_row)
  }

  # generate output for further modifications
  if (out_list == TRUE) {
    outlist <- list(vip = vip_plot, pdp = pdp_curves)
    return(outlist)
  }
}
