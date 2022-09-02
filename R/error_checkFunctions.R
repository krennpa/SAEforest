# Throwing general class error for generic functions -------------------------------------

class_error <- function(object) {
  if (!inherits(object, "SAEforest")) {
    stop("The object has to be created by the SAEforest package for 'SAEforest'- methods to work.")
  }
}


# Checking inputs for the main function SAEforest_model ----------------------------------

input_checks_model <- function(Y, X, dName, smp_data, pop_data, MSE, meanOnly, aggData, smearing,
                               popnsize, importance, OOsample_obs, ADDsamp_obs, w_min, B, B_adj,
                               B_MC, threshold, custom_indicator, initialRandomEffects, ErrorTolerance,
                               MaxIterations, na.rm) {
  if (!is.numeric(Y) || !data.frame(Y) %in% smp_data) {
    stop("Y must be a continuous vector containing the target variable. Additionally Y must be included in the data frame of survey sample data. See also help(SAEforest_model)")
  }

  if (sum(!(X %in% smp_data)) != 0) {
    stop("X specifies the explanatory variabels from the sample data set and must be included in the survey sample set. See also help(SAEforest_model)")
  }

  if (!is.data.frame(smp_data) || !is.data.frame(pop_data)) {
    stop("smp_data must be a data frame containing the survey sample data. pop_data must be a data frame containing the population data. See also help(SAEforest_model).")
  }

  if (!is.character(dName) || length(dName) != 1) {
    stop("dName must be an input of class character determining the variable name of the domain of interest provided in the survey as well as the population data. See also help(SAEforest_model).")
  }

  if (data.frame(smp_data[[dName]]) %in% X) {
    stop("X must not contain the domain-level identifier. X must contain only covariate data used for the estimation of fixed effects part of the model.")
  }

  if ((aggData == FALSE) && dim(pop_data)[1] < dim(smp_data)[1]) {
    stop("The population data set cannot have less observations than the sample data set.")
  }

  if (is.null(smp_data[[dName]]) || is.null(pop_data[[dName]])) {
    stop(paste('The survey sample data and the population data must contain information on domains. Both data frames must contain a column labelled by the same name specified under the input of "dName"'))
  }

  if (!(inherits(na.rm, "logical") && length(na.rm) == 1) || na.rm == FALSE && (sum(is.na(smp_data)) != 0 ||
    sum(is.na(pop_data)) != 0)) {
    stop("The survey data or the population data contain missing values. Please set na.rm = TRUE or consider imputation methods handling missing values in your data.")
  }

  if (!all(unique(as.character(na.omit(smp_data[[dName]]))) %in%
    unique(as.character(pop_data[[dName]])))) {
    stop("The survey sample data contains domains that are not contained in the population data.")
  }

  if (is.null(colnames(smp_data)) || is.null(colnames(pop_data)) ||
    sum((colnames(smp_data) %in% colnames(pop_data))) == 0) {
    stop("smp_data and pop_data must contain columnames for covariates that allow for a clear linkage of covariate data from the survey sample and the population level information. See also help(SAEforest_model)")
  }

  if (!sum(colnames(X) %in% colnames(pop_data)) == dim(X)[2]) {
    stop("The survey sample covariates X contain predictors that are not contained in the population data. Please ensure a clear linkage of covariate data from the survey sample and the population level information.")
  }

  if (!is.numeric(initialRandomEffects) || (length(initialRandomEffects) != 1 && length(initialRandomEffects) != length(Y))) {
    stop(paste("initialRandomEffects specify initional values of random effects for the MERF. Acceptable inputs are single values such as the default of 0 or numeric vectors of length: ", length(Y)))
  }

  if (!is.numeric(MaxIterations) || length(MaxIterations) != 1 || MaxIterations < 2) {
    stop("MaxIterations needs to be a single integer value, determining the number of maximum iterations for the convergence of the MERF algorithm. The value must be at least 2. See also help(SAEforest_model).")
  }

  if (!is.numeric(ErrorTolerance) || length(ErrorTolerance) != 1 || ErrorTolerance <= 0) {
    stop("ErrorTolerance needs to be a single integer value, determining the tolerance to monitor the convergence of the MERF algorithm. The value must be greater than 0. See also help(SAEforest_model).")
  }

  if (is.null(MSE) || !(MSE == "none" || MSE == "wild" || MSE == "nonparametric")) {
    stop('The options for MSE are "none" and "nonparametric". For nonlinear indicators (meanOnly = FALSE) the additional option "wild" exists. See also help(SAEforest_model).')
  }

  if (MSE != "none" && !(is.numeric(B) && length(B) == 1 && B > 1)) {
    stop("If MSE estimation is specified, B needs to be a single integer value, determining the number of MSE-bootstrap replications. The value must be larger than 1. See also help(SAEforest_model).")
  }

  if (!(inherits(aggData, "logical") && length(aggData) == 1)) {
    stop("The option aggData is logical. Insert TRUE if you only have aggregated covariate population data. See also help(SAEforest_model)")
  }

  if (!(inherits(meanOnly, "logical") && length(meanOnly) == 1)) {
    stop("The option meanOnly is logical. Insert TRUE if you are only interested in the estimation of domain-level means to save computation time. See also help(SAEforest_model)")
  }

  if (is.null(importance) || !(importance == "none" || importance == "impurity" || importance == "impurity_corrected" || importance == "permutation")) {
    stop('Variable importance is needed for vip plots. To reduce runtime, importance can be set to "none". Furhter options are "impurity", "impurity_corrected" or "permutation". See details on the variable importance mode with help(ranger).')
  }

  if ((aggData == TRUE) && (is.null(importance) || !(importance == "impurity" || importance == "impurity_corrected" || importance == "permutation"))) {
    stop('The calculation of calibration weights requires a concept of variable importance. The options are passed to the MERF and must be "impurity", "impurity_corrected" or "permutation". See details on the variable importance mode with help(ranger). See details on the caluclation of calibration weights help(SAEforest_model).')
  }

  if (!(inherits(smearing, "logical") && length(smearing) == 1)) {
    stop("The option smearing is logical. Insert TRUE for smearing based estimation and the value FALSE for MC-based estimation. See also help(SAEforest_model)")
  }

  if (MSE != "none" && !(is.numeric(B_adj) && length(B_adj) == 1 && B_adj > 1)) {
    stop("If MSE estimation is specified, B_adj needs to be a single integer value, determining the number of bootstrap replications for the adjustemt of residual variance. The value must be larger than 1. See also help(SAEforest_model).")
  }

  if (smearing != TRUE && !(is.numeric(B_MC) && length(B_MC) == 1 && B_MC > 1)) {
    stop("If MC-based estimation is specified, B_MC needs to be a single integer value, determining the number of generated bootstrap populations. The value must be larger than 1. See also help(SAEforest_model).")
  }

  if (!is.null(threshold) && !(is.numeric(threshold) && length(threshold) == 1) && !inherits(threshold, "function")) {
    stop("threshold needs to be a single numeric value or a function of Y. If it is NULL 60% of the median of Y is selected as threshold. See also help(SAEforest_model).")
  }

  if (inherits(threshold, "function") && !all(attributes(formals(threshold))$names == c("Y"))) {
    stop("If threshold is a function, the argument must depend exclusively on Y. Also, a single numeric value is possible as threshold. If threshld = NULL, 60% of the median of Y is selected as threshold. See also help(SAEforest_model).")
  }

  if (!is.null(custom_indicator)) {
    if (!inherits(custom_indicator, "list")) {
      stop("Additional indicators must be added via custom_indicator as a list of functions. For help see the example in help(SAEforest_model).")
    }

    N_custom <- length(custom_indicator)
    for (i in seq_len(N_custom)) {
      if (!inherits(custom_indicator[[i]], "function")) {
        stop("All elements of the list need to be functions. These functions for custom indicators must contain exactly two arguments: Y and threshold; threshold might be NULL. For help see the Example in help(SAEforest_model).")
      } else if (inherits(custom_indicator[[i]], "function") &&
        !all(names(formals(custom_indicator[[i]])) == c("Y", "threshold"))) {
        stop("Functions for custom indicators must contain exactly two arguments: Y and threshold; threshold might be NULL. For help see the example in help(SAEforest_model).")
      }
    }
  }

  if ((aggData == TRUE) && (MSE != "none" && (!is.data.frame(popnsize) || is.null(popnsize[[dName]]) || dim(popnsize)[1] != dim(pop_data)[1] ||
    dim(popnsize)[2] > 2 || !is.numeric(popnsize[, !colnames(popnsize) %in% dName]) || !all.equal(as.character(popnsize[[dName]]), as.character(pop_data[[dName]]))))) {
    stop(paste("popnsize must be a data frame with two columns: One column named ", dName, "must contain the domains specifiers. The other column must contain information on the population size of domains. See also help(SAEforest_model)."))
  }

  if ((aggData == TRUE) && (!is.numeric(OOsample_obs) || length(OOsample_obs) != 1 || OOsample_obs < 0)) {
    stop('OOsample_obs needs to be a single integer value, determining the amount of observations sampled from the "closest" area for out-of-sample areas. See also help(SAEforest_model).')
  }

  if ((aggData == TRUE) && (!is.numeric(ADDsamp_obs) || length(ADDsamp_obs) != 1 || ADDsamp_obs < 0)) {
    stop('ADDsamp_obs needs to be a single integer value, determining the amount of observations sampled from the "closest" area in the case of failure of calculation of calibration weights. See also help(SAEforest_model).')
  }

  if ((aggData == TRUE) && (!is.numeric(w_min) || length(w_min) != 1 || w_min < 2 || w_min > dim(pop_data)[2])) {
    stop("w_min needs to be a single integer value, determining the minimum amount of covariates incorporating auxilliary information for the assessment of calibration weights. Thus, w_min must be smaller or equal to the number of existing covariates. See also help(SAEforest_model).")
  }

  if ((aggData == TRUE || meanOnly == TRUE) && (is.null(MSE) || !(MSE == "none" || MSE == "nonparametric"))) {
    stop('The options for the MSE for means are "none" or "nonparametric".')
  }
}

# Checking inputs for the generic plot function ------------------------------------------

input_checks_plot <- function(num_features, alpha, include_type, horizontal,
                              lsize, grid_row, out_list, pdp_plot, gg_theme) {
  if (!is.numeric(num_features) || length(num_features) != 1 || num_features < 2) {
    stop("num_features needs to be a single integer value, determining the amount of covariates for which variable importance should be plotted.")
  }

  if (!is.numeric(alpha) || length(alpha) != 1) {
    stop("alpha needs to be a single integer value between 0 and 1, determining the transperancy of colors.")
  }

  if (!is.numeric(grid_row) || length(grid_row) != 1) {
    stop("grid_row must be a single integer value determining the number of rows for the joint plot of pdps.")
  }

  if (!is.logical(include_type)) {
    stop("include_type must be logical! See details help(plot.SAEforest).")
  }

  if (!is.logical(horizontal)) {
    stop("horizontal must be logical! See details help(plot.SAEforest).")
  }

  if (!is.logical(out_list)) {
    stop("out_list must be logical! See details help(plot.SAEforest).")
  }

  if (!is.logical(pdp_plot)) {
    stop("pdp_plot must be logical! See details help(plot.SAEforest).")
  }

  if (!is.numeric(lsize) || length(lsize) != 1) {
    stop("lsize must be a single integer value determining the line size for pdp plots.")
  }
  if (!inherits(gg_theme, "gg") || length(lsize) != 1) {
    stop("gg_theme must be a single function of type 'gg' determining the theme of plots.")
  }
}


# Checking inputs for summarize_indicators function --------------------------------------

summarize_indicators_check <- function(object, indicator, MSE, CV) {
  class_error(object)

  if (is.null(object$MSE_Estimates) && (MSE == TRUE || CV == TRUE)) {
    stop("No MSE estimates in SAEforest object: arguments MSE and CV have to be FALSE or a new SAEforest object with variance/MSE needs to be generated.")
  }
  if (!(inherits(MSE, "logical") && length(MSE) == 1)) {
    stop("MSE must be a logical value. Set MSE to TRUE or FALSE.")
  }
  if (!(inherits(CV, "logical") && length(CV) == 1)) {
    stop("CV must be a logical value. Set CV to TRUE or FALSE.")
  }
  if ((inherits(object, "SAEforest_mean") || inherits(object, "SAEforest_meanAGG")) && !inherits(indicator, "Mean") && (!inherits(indicator, "all") ||
     !inherits(indicator, "All"))) {
    warning('For objects of class SAEforest_mean or SAEforest_meanAGG, only results for the Indicator "Mean" are reported')
  }
  if (inherits(object, "SAEforest_nonLin")) {
    if (is.null(indicator) || !all(indicator == "all" | indicator == "All" | indicator == "Mean" |
      indicator == "Quant10" | indicator == "Quant25" | indicator == "Median" |
      indicator == "Quant75" | indicator == "Quant90" | indicator == "Gini" |
      indicator == "Hcr" | indicator == "Pgap" | indicator == "Qsr" |
      indicator %in% names(object$Indicators[-1]))) {
      stop(paste0("The argument indicator is set to ", indicator, ". The argument only allows to be set to all, a name of estimated indicators or indicator groups defined by c(). For details see help(summarize_indicators)"))
    }
  }
}

# Checking inputs for map_indicators function --------------------------------------------

map_indicators_check <- function(object, indicator, MSE, CV, map_obj, map_dom_id, map_tab, color,
                                 return_data, return_plot, gg_theme) {
  class_error(object)

  if (!(inherits(CV, "logical") && length(CV) == 1)) {
    stop("CV must be a logical value. Set CV to TRUE or FALSE.")
  }
  if (!(inherits(MSE, "logical") && length(MSE) == 1)) {
    stop("MSE must be a logical value. Set MSE to TRUE or FALSE.")
  }
  if (!(inherits(return_data, "logical") && length(return_data) == 1)) {
    stop("return_data must be a logical value. Set return_data to TRUE or FALSE.")
  }
  if (!(inherits(return_plot, "logical") && length(return_plot) == 1)) {
    stop("return_plot must be a logical value. Set return_plot to TRUE or FALSE.")
  }
  if (is.null(object$MSE_Estimates) && (MSE == TRUE || CV == TRUE)) {
    stop("No MSE estimates in SAEforest object: arguments MSE and CV have to be FALSE or a new SAEforest object with variance/MSE needs to be generated.")
  }
  if (length(color) != 2 || !is.vector(color)) {
    stop("col needs to be a vector of length 2 defining the starting and upper color of the map-plot")
  }
  if (!(map_dom_id %in% names(map_obj))) {
    stop(paste0(map_dom_id, " is not contained in map_obj. Please provide valid variable name for pop_domains."))
  }

  if (inherits(object, "SAEforest_nonLin")) {
    if (is.null(indicator) || !all(indicator == "all" | indicator == "All" | indicator == "Mean" |
      indicator == "Quant10" | indicator == "Quant25" | indicator == "Median" |
      indicator == "Quant75" | indicator == "Quant90" | indicator == "Gini" |
      indicator == "Hcr" | indicator == "Pgap" | indicator == "Qsr" |
      indicator %in% names(object$Indicators[-1]))) {
      stop(paste0("The argument indicator is set to ", indicator, ". The argument only allows to be set to all, a name of estimated indicators or indicator groups defined by c(). For details see help(summarize_indicators)"))
    }
  }
  if (is.null(object$MSE_Estimates) && (MSE == TRUE || CV == TRUE)) {
    stop("No MSE estimates in SAEforest object: arguments MSE and CV have to be FALSE or a new SAEforest object with variance/MSE needs to be generated.")
  }

  if (!inherits(gg_theme, "gg")) {
    stop("gg_theme must be a single function of type 'gg' determining the theme of plots.")
  }

  if (!inherits(map_obj, "SpatialPolygonsDataFrame") ||
    attr(class(map_obj), "package") != "sp") {
    stop("map_obj is not of class SpatialPolygonsDataFrame from the sp package")
  }
  if (!is.null(map_tab) && !(inherits(map_tab, "data.frame") &&
    dim(map_tab)[2] == 2)) {
    stop("If the domain-identifiers in the data object and shape file differ, a mapping table needs to be used. This table needs to be a data frame with two colums. See also help(map_indicators).")
  }
}

# Checking inputs for MERFranger function ------------------------------------------------

input_checks_MERF <- function(Y, X, data, initialRandomEffects,
                              ErrorTolerance, MaxIterations, importance, na.rm) {
  if (!is.numeric(Y) || !length(Y) == dim(data)[1]) {
    stop("Y must be a metric vector containing the target variable. Additionally Y must be included in the data frame of survey sample data. See also help(MERFranger)")
  }

  if (sum(!(X %in% data)) != 0) {
    stop("X specifies the explanatory variabels from the sample data set and must be included in the survey sample set. See also help(MERFranger)")
  }

  if (!is.data.frame(data)) {
    stop("smp_data must be a data frame containing the survey sample data. See also help(MERFranger).")
  }

  if (!(inherits(na.rm, "logical") && length(na.rm) == 1) || na.rm == FALSE && (sum(is.na(data)) != 0)) {
    stop("The survey data contains missing values. Please set na.rm = TRUE or consider imputation methods handling missing values in your data.")
  }

  if (!is.numeric(initialRandomEffects) || (length(initialRandomEffects) != 1 && length(initialRandomEffects) != length(Y))) {
    stop(paste("initialRandomEffects specify initional values of random effects for the MERF. Acceptable inputs are single values such as the default of 0 or numeric vectors of length: ", length(Y)))
  }

  if (!is.numeric(MaxIterations) || length(MaxIterations) != 1 || MaxIterations < 2) {
    stop("MaxIterations needs to be a single integer value, determining the number of maximum iterations for the convergence of the MERF algorithm. The value must be at least 2. See also help(MERFranger).")
  }

  if (!is.numeric(ErrorTolerance) || length(ErrorTolerance) != 1 || ErrorTolerance <= 0) {
    stop("ErrorTolerance needs to be a single integer value, determining the tolerance to monitor the convergence of the MERF algorithm. The value must be greater than 0. See also help(MERFranger).")
  }

  if (is.null(importance) || !(importance == "none" || importance == "impurity" || importance == "impurity_corrected" || importance == "permutation")) {
    stop('Variable importance is needed for vip plots. To reduce runtime, importance can be set to "none". Furhter options are "impurity", "impurity_corrected" or "permutation". See details on the variable importance mode with help(ranger).')
  }
}


# Checking inputs for tune_parameters function -------------------------------------------

input_checks_tune <- function(Y, X, data, seed,
                              gg_theme, plot_res, return_plot) {
  if (!is.numeric(Y) || !length(Y) == dim(data)[1]) {
    stop("Y must be a metric vector containing the target variable. Additionally Y must be included in the data frame of survey sample data. See also help(tune_parameters)")
  }

  if (sum(!(X %in% data)) != 0) {
    stop("X specifies the explanatory variabels from the sample data set and must be included in the survey sample set. See also help(tune_parameters)")
  }

  if (!is.data.frame(data)) {
    stop("smp_data must be a data frame containing the survey sample data. See also help(tune_parameters).")
  }

  if (!inherits(gg_theme, "gg")) {
    stop("gg_theme must be a single function of type 'gg' determining the theme of plots.")
  }
  if (!(inherits(return_plot, "logical") && length(return_plot) == 1)) {
    stop("return_plot must be a logical value. Set return_plot to TRUE or FALSE.")
  }
  if (!(inherits(plot_res, "logical") && length(plot_res) == 1)) {
    stop("plot_res must be a logical value. Set plot_res to TRUE or FALSE.")
  }
  if (!is.numeric(seed) || length(seed) != 1) {
    stop("seed needs to be a single integer value, determining a seed to be set for reproduceability.")
  }
}
