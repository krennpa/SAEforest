#' Simulated EU-SILC data - aggregated population data
#'
#' The data set includes synthetic EU-SILC data and is taken from the package \pkg{emdi}.
#' Originally, the data builds on \code{eusilcP} from package \pkg{simFrame}
#' and is reduced to 15 variables including district identifiers as well as aggregated
#' household level covariates. Therefore, except for the variables \code{ratio_n} and
#' \code{district}, the variables are mean values per district.
#'
#'
#' @format A data frame with 94 observations and 15 variables:
#' \describe{
#' \item{eqsize}{numeric; the equivalized household size according to the
#' modified OECD scale.}
#' \item{cash}{numeric; employee cash or near cash income (net).}
#' \item{self_empl}{numeric; cash benefits or losses from self-employment (net).}
#' \item{unempl_ben}{numeric; unemployment benefits (net).}
#' \item{age_ben}{numeric; old-age benefits (net).}
#' \item{surv_ben}{numeric; survivor's benefits (net).}
#' \item{sick_ben}{numeric; sickness benefits (net).}
#' \item{dis_ben}{numeric; disability benefits (net).}
#' \item{rent}{numeric; income from rental of a property or land (net).}
#' \item{fam_allow}{numeric; family/children related allowances (net).}
#' \item{house_allow}{numeric; housing allowances (net).}
#' \item{cap_inv}{numeric; interest, dividends, profit from capital investments
#'  in unincorporated business (net).}
#' \item{tax_adj}{numeric; repayments/receipts for tax adjustment (net).}
#' \item{ratio_n}{numeric; ratios of the population size per area and the total population size.}
#' \item{district}{factor; Austrian districts (94 levels).}
#' }
#' @docType data
"eusilcA_popAgg"
