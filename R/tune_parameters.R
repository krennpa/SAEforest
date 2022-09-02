#' Tuning and cross-validation of MERF parameters
#'
#' Function \code{tune_parameters} allows to tune parameters for the implemented MERF method. Essentially,
#' this function can be understood as a modified wrapper for \link[caret]{train} from the package \pkg{caret},
#' treating MERFs as a custom method.
#'
#' @details Tuning can be performed on the following four parameters: \code{num.trees} (the number of trees
#' for a forest), \code{mtry} (number of variables as split candidates at in each node), \code{min.node.size}
#' (minimal individual node size) and \code{splitrule} (general splitting rule). For details see
#' \link[ranger]{ranger}.
#'
#' @param Y Continuous input value of target variable.
#' @param X Matrix or data.frame of predictive covariates.
#' @param data data.frame of survey sample data including the specified elements of \code{Y} and
#' \code{X}.
#' @param dName Character specifying the name of domain identifier, for which random intercepts
#' are modeled.
#' @param trControl Control parameters passed to \link[caret]{train}. Most important
#' parameters are \code{method} ("repeatedcv" for x-fold cross-validation), \code{number} (the number of folds)
#' and \code{repeats} (the number of repetitions). For further details see \link[caret]{trainControl}
#' and the example below.
#' @param tuneGrid A data.frame with possible tuning values. The columns must have the same names as the
#' tuning parameters. For this tuning function the grid must comprise entries for the following parameters:
#' \code{num.trees, mtry, min.node.size, splitrule}.
#' @param seed Enabling reproducibility of for cross-validation and tuning. Defaults to \code{11235}.
#' @param gg_theme Specify a predefined theme from \pkg{ggplot2}. Defaults to \code{theme_minimal}.
#' @param plot_res Optional logical. If \code{TRUE}, the plot with results of cross-validation and tuning
#' is shown. Defaults to \code{TRUE}.
#' @param return_plot If set to \code{TRUE}, a list of the comparative plot produced by \pkg{ggplot2}
#' is returned for further individual customization and processing.
#' @param na.rm Logical. Whether missing values should be removed. Defaults to \code{TRUE}.
#' @param ... Additional parameters are directly passed to the random forest \link[ranger]{ranger} and/or
#' the training function \link[caret]{train}. For further details on possible parameters and examples
#' see \link[ranger]{ranger} or \link[caret]{train}.
#'
#' @return Prints requested optimal tuning parameters and (if requested) an additional
#' comparative plot produced by \pkg{ggplot2}.
#'
#' @seealso \code{\link{SAEforest}}, \code{\link{MERFranger}}, \code{\link[caret]{train}},
#' \code{\link[ggplot2]{ggplot}}
#'
#' @examples
#' \donttest{
#' # Loading data
#' data("eusilcA_pop")
#' data("eusilcA_smp")
#' library(caret)
#'
#' income <- eusilcA_smp$eqIncome
#' X_covar <- eusilcA_smp[, -c(1, 16, 17, 18)]
#'
#' # Specific characteristics of Cross-validation
#' fitControl <- trainControl(method = "repeatedcv", number = 5,
#'                            repeats = 1)
#'
#' # Define a tuning-grid
#' merfGrid <- expand.grid(num.trees = 50, mtry = c(3, 7, 9),
#'                         min.node.size = 10, splitrule = "variance")
#'
#' tune_parameters(Y = income, X = X_covar, data = eusilcA_smp,
#'                 dName = "district", trControl = fitControl,
#'                 tuneGrid = merfGrid)
#' }
#'
#' @export
#' @import caret
#' @importFrom ggplot2 theme_minimal ggplot

tune_parameters <- function(Y,
                            X,
                            data,
                            dName,
                            trControl,
                            tuneGrid,
                            seed = 11235,
                            gg_theme = theme_minimal(),
                            plot_res = TRUE,
                            return_plot = FALSE,
                            na.rm = TRUE,
                            ...) {
  input_checks_tune(
    Y = Y, X = X, data = data, seed = seed, gg_theme = gg_theme,
    plot_res = plot_res, return_plot = return_plot
  )

  if (na.rm == TRUE) {
    comp_smp <- complete.cases(data)
    data <- data[comp_smp, ]
    Y <- Y[comp_smp]
    X <- X[comp_smp, ]
  }

  MERF <- define_method(X = X, dName = dName, ...)

  set.seed(seed)
  MERFtune <- caret::train(
    x = cbind(X, data[dName]),
    y = Y,
    method = MERF,
    trControl = trControl,
    tuneGrid = tuneGrid, ...
  )
  print(MERFtune)

  paramPlot <- ggplot(MERFtune) + gg_theme

  if (plot_res) {
    print(paramPlot)
  }
  if (return_plot) {
    return(paramPlot)
  }
}


# Defining a custom MERF method for caret framework ---------------------------------------
define_method <- function(X, dName, ...) {

  random <- paste0(paste0("(1|", dName), ")")
  colnames_in <- c("Y", colnames(X), dName)
  m_example <- c(sqrt(dim(X)[2]), dim(X)[2])

  MERF <- list(
    type = "Regression",
    library = NULL,
    loop = NULL
  )
  # parameters to be tuned
  tp <- data.frame(
    parameter = c("num.trees", "mtry", "min.node.size", "splitrule"),
    class = c("numeric", "numeric", "numeric", "character"),
    label = c("num.trees", "mtry", "min.node.size", "splitrule")
  )

  MERF$parameters <- tp

  # default grid
  MERFGrid <- function(x, y, len = NULL, search = "grid") {
    if (search == "grid") {
      out <- expand.grid(
        num.trees = c(500, 1000),
        mtry = m_example,
        min.node.size = 5, splitrule = "variance"
      )
    } else {
      stop("random search not yet implemented")
    }
    out
  }
  # append to list
  MERF$grid <- MERFGrid

  MERFFit <- function(x, y, wts, param, lev, last, weights, classProbs, ...) {
    dat <- as.data.frame(cbind(y, x))
    colnames(dat) <- colnames_in

    MERFranger(
      X = x[, (2:dim(dat)[2] - 1)],
      Y = dat[, 1],
      num.trees = param$num.trees,
      mtry = param$mtry,
      min.node.size = param$min.node.size,
      splitrule = param$splitrule,
      random = random, data = dat,
      ...
    )
  }
  # append to list
  MERF$fit <- MERFFit

  MERFPred <- function(modelFit, newdata, preProc = NULL, submodels = NULL) {
    predict.MERFranger(modelFit, newdata)
  }

  MERF$predict <- MERFPred
  MERF$prob <- 0
  MERF$Sort <- function(x) {
    x
  }

  return(MERF)
}
