#' Visualizes regional disaggregated estimates on a map
#'
#' Function \code{map_indicators} visualizes estimates from a
#' \code{\link{SAEforestObject}} on a specified map. The function can be seen as a modified
#' wrapper of \code{\link[emdi]{map_plot}} from the package \pkg{emdi}.
#'
#' @param object An object of class \code{SAEforest}, containing estimates to be visualized.
#' @param indicator Optional character vector specifying indicators to be mapped: (i)
#' all calculated indicators ("all"); (ii) each default indicators name: "Mean",
#' "Quant10", "Quant25", "Median", "Quant75", "Quant90", "Gini", "Hcr", "Pgap", "Qsr"
#' or the function name/s of "custom_indicator/s"; (iii) a vector of names of indicators.
#' If the \code{object} is estimated by \code{\link{SAEforest_mean}},
#' indicator arguments are ignored and only "Mean" is visualized.
#' @param MSE Logical. If \code{TRUE}, the MSE is also visualized.
#' Defaults to \code{FALSE}.
#' @param CV Logical. If \code{TRUE}, the CV is also visualized.
#' Defaults to \code{FALSE}.
#' @param map_obj An \code{SpatialPolygonsDataFrame} object as defined by the
#' \pkg{sp} package on which the data should be visualized.
#' @param map_dom_id Character string containing the name of a variable in
#' \code{map_obj} that indicates the domains.
#' @param map_tab A \code{data.frame} object with two columns that matches the
#' domain variable from the population data set (first column) with the domain
#' variable in the map_obj (second column). This should only be used if domain level identifiers
#' are different in both objects.
#' @param color A \code{vector} of length 2 defining the lowest and highest color in the map.
#' @param scale_points A structure defining the lowest, the mid and the highest
#' value of the colorscale. If a numeric vector of length two is given, this scale
#' will be used for every plot. Alternatively, a list defining colors for each plot separately may be given.
#' @param guide Character passed to \code{scale_colour_gradient} from \pkg{ggplot2}.
#' Possible values are "none", "colourbar", and "legend".
#' @param return_data If set to \code{TRUE}, a fortified data frame including the
#' map data as well as the chosen indicators is returned. Customized maps can
#' easily be obtained from this data frame via the package \pkg{ggplot2}. Defaults to \code{FALSE}.
#' @param return_plot If set to \code{TRUE}, a list of individual plots produced by \pkg{ggplot2}
#' is returned for further individual customization and processing
#' @param gg_theme Specify a predefined theme from \pkg{ggplot2}. Default is set to \code{theme_minimal}.
#'
#' @return Creates required plots and if selected, a fortified data.frame and a list of plots.
#'
#' @seealso \code{\link{SAEforest}}, \code{\link[maptools]{readShapePoly}},
#' \code{\link[sp]{SpatialPolygonsDataFrame}}, \code{\link[ggplot2]{ggplot}}.
#'
#' @examples
#' \dontrun{
#'#Loading data
#'data("eusilcA_pop")
#'data("eusilcA_smp")
#'
#'income <- eusilcA_smp$eqIncome
#'X_covar <- eusilcA_smp[,-c(1,16,17,18)]
#'
#'#Example 1:
#'#Calculating point estimates and discussing basic generic functions
#'
#'model1 <- SAEforest_mean(Y = income, X = X_covar, dName = "district",
#'                        smp_data = eusilcA_smp, pop_data = eusilcA_pop, num.trees=50)
#'
#'# Load shape file
#'load_shapeaustria()
#'
#'# Create map plot for mean indicator - point and MSE estimates but no CV
#'map_indicators(object = model1, MSE = FALSE, CV = FALSE,
#'          map_obj = shape_austria_dis, indicator = c("Mean"),
#'          map_dom_id = "PB")
#'
#'# Create a suitable mapping table to use numerical identifiers of the shape
#'# file
#'
#'# First find the right order
#'dom_ord <- match(shape_austria_dis@data$PB, model1$Indicators$district)
#'
#'# Create the mapping table based on the order obtained above
#'map_tab <- data.frame(pop_data_id = model1$Indicators$district[dom_ord],
#'                     shape_id = shape_austria_dis@data$BKZ)
#'
#'# Create map plot for mean indicator - using the numerical domain
#'identifiers of the shape file. Additionally save the figure in as a list element.
#'
#'map_obj <- map_indicators(object = model1, MSE = FALSE, CV = FALSE,
#'            map_obj = shape_austria_dis, indicator = c("Mean"),
#'           map_dom_id = "BKZ", map_tab = map_tab, return_plot = TRUE)
#'}
#'
#' @export
#' @importFrom reshape2 melt
#' @importFrom ggplot2 aes geom_polygon facet_wrap fortify coord_equal labs
#' @importFrom ggplot2 theme element_blank guides scale_fill_gradient
#' @importFrom ggplot2 scale_colour_gradient theme_minimal ggplot
#' @import maptools

map_indicators <- function(object,
                     indicator = "all",
                     MSE = FALSE,
                     CV = FALSE,
                     map_obj = NULL,
                     map_dom_id = NULL,
                     map_tab = NULL,
                     color = c("white", "darkgreen"),
                     scale_points = NULL,
                     guide = "colourbar",
                     return_data = FALSE,
                     return_plot = FALSE,
                     gg_theme = theme_minimal()){

  map_indicators_check(object = object, indicator = indicator, MSE = MSE, CV = CV, map_obj = map_obj,
                       map_dom_id = map_dom_id, map_tab = map_tab, color = color,
                       return_data = return_data, return_plot = return_plot, gg_theme = gg_theme)


  if (is.null(map_obj)) {
    message("No Map Object has been provided. An artificial polygone is used for
            visualization")
    map_pseudo(object = object , indicator = indicator, panelplot = FALSE, MSE = MSE,
               CV = CV, gg_theme = gg_theme, return_plot = return_plot)
  }

    plot_real(object,
              indicator = indicator,
              MSE = MSE,
              CV = CV,
              map_obj = map_obj,
              map_dom_id = map_dom_id,
              map_tab = map_tab,
              col = color,
              scale_points = scale_points,
              return_data = return_data,
              guide = guide, gg_theme = gg_theme, return_plot = return_plot
    )
  }

map_pseudo <- function(object, indicator, panelplot, MSE, CV, gg_theme, return_plot)
{
  x <- y <- id <- value <- NULL #avoid note due to usage in ggplot
  values <-  summarize_indicators(object = object, indicator = indicator,
                        MSE = MSE, CV = CV)$ind
  indicator <- colnames(values)[-1]

  tplot <- get_polygone(values = values)

  plot_list <- vector(mode="list", length = length(indicator))
  names(plot_list) <- indicator

  if (panelplot) {
    ggplot(tplot, aes(x = x, y = y)) + geom_polygon(aes(group=id, fill = value))
    + facet_wrap( ~ variable, ncol = ceiling(sqrt(length(unique(tplot$variable))))) + gg_theme
  } else {
    for (ind in indicator) {
      plot_list[[ind]] <- eval(substitute(ggplot(tplot[tplot$variable == ind,], aes(x = x, y = y)) +
        ggtitle(paste0(ind)) + geom_polygon(aes(group = id, fill = value)) + gg_theme, list(ind=ind)))

      print(plot_list[[ind]])
      cat("Press [enter] to continue")
      line <- readline()
    }
  }
  if(return_plot){
    return(plot_list)
  }
}

plot_real <- function(object,
                      indicator = "all",
                      MSE = FALSE,
                      CV = FALSE,
                      map_obj = NULL,
                      map_dom_id = NULL,
                      map_tab = NULL,
                      col = col,
                      scale_points = NULL,
                      return_data = FALSE,
                      guide = NULL,
                      gg_theme,
                      return_plot =FALSE
) {
  if (!is.null(map_obj) && is.null(map_dom_id)) {
    stop("No correct ID for the map object is given")
  }
  long <- lat <- group <- NULL

  map_data <- summarize_indicators(object = object, indicator = indicator,
                         MSE = MSE, CV = CV)$ind

  if (!is.null(map_tab)) {
    map_data <- merge(x = map_data, y = map_tab,
                      by.x = names(object$Indicators)[1], by.y = names(map_tab)[1])
    matcher <- match(map_obj@data[map_dom_id][,1], map_data[,names(map_tab)[2]])

    if (any(is.na(matcher))) {
      if (all(is.na(matcher))) {
        stop("Domain identifier of map_tab and Map object do not match. Check map_tab")
      } else {
        warnings("Not all domain identifiers of map_tab and Map object could be matched.
                 Check map_tab")
      }
    }
    map_data <- map_data[matcher,]
    map_data <- map_data[,!colnames(map_data) %in% c(names(object$Indicators)[1],
                                                     map_dom_id,
                                                     names(map_tab)), drop = F]
  } else {
    matcher <- match(map_obj@data[map_dom_id][,1], map_data[,names(object$Indicators)[1]])

    if (any(is.na(matcher))) {
      if (all(is.na(matcher))) {
        stop("dName of SAEforest object and Map object do not match. Try using map_tab")
      } else {
        warnings("Not all domain identifiers of EMDI and Map object could be matched.
                 Try using map_tab")
      }
    }
    map_data <- map_data[matcher, ]
  }

  map_obj@data[colnames(map_data)] <- map_data


  map_obj.fort <- fortify(map_obj, region = map_dom_id)
  map_obj.fort <- merge(map_obj.fort, map_obj@data,
                        by.x = "id", by.y = map_dom_id)

  indicator <- colnames(map_data)
  indicator <- indicator[!(indicator %in% names(object$Indicators)[1])]

  plot_list <- vector(mode="list", length = length(indicator))
  names(plot_list) <- indicator

  for (ind in indicator) {
    map_obj.fort[ind][,1][!is.finite(map_obj.fort[ind][,1])] <- NA
    scale_point <- get_scale_points(map_obj.fort[ind][,1], ind, scale_points)
    plot_list[[ind]] <- eval(substitute(ggplot(map_obj.fort, aes(long, lat, group = group,
                                                 fill = map_obj.fort[ind][,1])) +
      geom_polygon(color = "azure3") + coord_equal() +
      labs(x = "", y = "", fill = ind) +
      ggtitle(gsub(pattern = "_",replacement = " ",x = ind)) +
      scale_fill_gradient(low = col[1], high = col[2],limits = scale_point,
                          guide = guide) + gg_theme +
      theme(axis.ticks = element_blank(), axis.text = element_blank(),
            legend.title = element_blank()), list(ind=ind)))

    print(plot_list[[ind]])

    if (!ind == tail(indicator,1)) {
      cat("Press [enter] to continue")
      line <- readline()
    }
  }
  if(return_data == TRUE && return_plot == FALSE) {
    return(map_obj.fort)
  }
  if(return_data == FALSE && return_plot == TRUE){
    return(plot_list)
  }
  if(return_data == TRUE && return_plot == TRUE){
    return(list(mapObj = map_obj.fort, plotObj = plot_list))
  }
}

# from emdi
get_polygone <- function(values) {
  if (is.null(dim(values))) {
    values = as.data.frame(values)
  }
  n <- nrow(values)
  cols <- ceiling(sqrt(n))
  n <- cols^2

  values["id"] <- seq_len(nrow(values))

  poly <- data.frame(id = rep(seq_len(n), each = 4),
                     ordering = seq_len((n*4)),
                     x = c(0,1,1,0) + rep(0:(cols - 1), each = (cols * 4)),
                     y = rep(c(0,0,1,1) + rep(0:(cols - 1), each = 4), cols)
  )

  combo <- merge(poly, values, by = "id", all = TRUE, sort = FALSE)
  reshape2::melt(combo[order(combo$ordering),], id.vars = c("id","x","y","ordering"))
}

# from emdi
get_scale_points <- function(y, ind, scale_points){
  result <- NULL
  if (!is.null(scale_points)) {
    if (class(scale_points) == "numeric" && length(scale_points) == 2) {
      result <- scale_points
    } else {
      splt <- strsplit(ind, "_\\s*(?=[^_]+$)", perl = TRUE)[[1]]
      indicator_name <- splt[1]
      if (length(splt) == 2) {
        measure <- splt[2]
      } else {
        measure <- "ind"
      }
      if (indicator_name %in% names(scale_points)) {
        pointset <- scale_points[[indicator_name]]
        try(result <- pointset[[measure]])
      }
      if (is.null(result) || length(result) != 2)
      {
        warning("scale_points is of no apropriate form, default values will
                 be used. See the descriptions and examples for details")
        result <- NULL
      }
    }
  }
  if (is.null(result)) {
    rg <- range(y, na.rm = TRUE)
    result <- rg
  }
  return(result)
}
