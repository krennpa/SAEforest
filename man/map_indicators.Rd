% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/map_indicators.R
\name{map_indicators}
\alias{map_indicators}
\title{Visualizes disaggregated estimates on a map}
\usage{
map_indicators(
  object,
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
  gg_theme = theme_minimal()
)
}
\arguments{
\item{object}{An object of class \code{SAEforest}, containing estimates to be visualized.}

\item{indicator}{Optional character vector specifying indicators to be mapped: (i)
all calculated indicators ("all"); (ii) default indicators name: "Mean",
"Quant10", "Quant25", "Median", "Quant75", "Quant90", "Gini", "Hcr", "Pgap", "Qsr"
or the function name/s of "custom_indicator/s"; (iii) a vector of names of indicators.
If the \code{object} is estimated with option \code{meanOnly = TRUE},
indicator arguments are ignored and only "Mean" is visualized.}

\item{MSE}{Logical. If \code{TRUE}, the MSE is also visualized.
Defaults to \code{FALSE}.}

\item{CV}{Logical. If \code{TRUE}, the CV is also visualized.
Defaults to \code{FALSE}.}

\item{map_obj}{An \code{"sf", "data.frame"} object as defined by the
\pkg{sf} package on which the data should be visualized.}

\item{map_dom_id}{Character string containing the name of a variable in
\code{map_obj} that indicates the domains.}

\item{map_tab}{A \code{data.frame} object with two columns that matches the
domain variable from the population data set (first column) with the domain
variable in the map_obj (second column). This should only be used if domain-level identifiers
are different in both objects.}

\item{color}{A \code{vector} of length 2 defining the lowest and highest color in the map.}

\item{scale_points}{A structure defining the lowest, the mid and the highest
value of the colorscale. If a numeric vector of length two is given, this scale
will be used for every plot. Alternatively, a list defining colors for each plot separately may be given.}

\item{guide}{Character passed to \code{scale_colour_gradient} from \pkg{ggplot2}.
Possible values are "none", "colourbar", and "legend".}

\item{return_data}{If set to \code{TRUE}, a fortified data frame including the
map data as well as the chosen indicators is returned. Customized maps can
easily be obtained from this data frame via the package \pkg{ggplot2}. Defaults to \code{FALSE}.}

\item{return_plot}{If set to \code{TRUE}, a list of individual plots produced by \pkg{ggplot2}
is returned for further individual customization and processing.}

\item{gg_theme}{Specify a predefined theme from \pkg{ggplot2}. Defaults to \code{theme_minimal}.}
}
\value{
Creates required plots and if selected, a fortified data.frame and a list of plots.
}
\description{
Function \code{map_indicators} visualizes estimates from a
\code{\link{SAEforestObject}} on a specified map. The function can be seen as a modified
wrapper of \code{map_plot} from the package \pkg{emdi}.
}
\examples{
\donttest{
# Loading data
data("eusilcA_pop")
data("eusilcA_smp")

# Load shape file
load_shapeaustria()

income <- eusilcA_smp$eqIncome
X_covar <- eusilcA_smp[, -c(1, 16, 17, 18)]

# Example 1:
# Calculating point estimates and discussing basic generic functions

model1 <- SAEforest_model(Y = income, X = X_covar, dName = "district",
                          smp_data = eusilcA_smp, pop_data = eusilcA_pop,
                          num.trees = 50)

# Create map plot for mean indicator - point and MSE estimates but no CV

map_indicators(object = model1, MSE = FALSE, CV = FALSE, map_obj = shape_austria_dis,
                          indicator = c("Mean"), map_dom_id = "PB")

# Create a suitable mapping table to use numerical identifiers of the shape
# file

# First find the right order
dom_ord <- match(shape_austria_dis$PB, model1$Indicators$district)

# Create the mapping table based on the order obtained above
map_tab <- data.frame(pop_data_id = model1$Indicators$district[dom_ord],
                      shape_id = shape_austria_dis$BKZ)

# Create map plot for mean indicator - using the numerical domain
# identifiers of the shape file. Additionally save the figure in as a list element.

map_obj <- map_indicators(object = model1, MSE = FALSE, CV = FALSE,
                          map_obj = shape_austria_dis, indicator = c("Mean"),
                          map_dom_id = "BKZ", map_tab = map_tab, return_plot = TRUE)
}

}
\seealso{
\code{\link{SAEforest}}, \code{\link[sf]{sf}}, \code{\link[ggplot2]{ggplot}}.
}
