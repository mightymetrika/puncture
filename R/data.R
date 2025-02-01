#' Power Analysis Results
#'
#' Results from power analysis with a 0.5 slope of interest and a 0.1 correlation
#' between predictors.
#'
#'Simulation framework modified from:
#'https://stats.stackexchange.com/questions/1866/how-to-simulate-a-custom-power-analysis-of-an-lm-model-using-r
#'Framework presented on stackexchange thread cites Bolker (2008) <ISBN: 9780691125220>
#'
#' @format ## `bb_pwer`
#' A data frame with 76 rows and 2 columns:
#' \describe{
#'   \item{N}{Sample size}
#'   \item{Pwr}{Rejection rate}
#' }
"bb_pwer"

#' Type 1 Error Results
#'
#' Results from power analysis with a 0 slope of interest and a 0.1 correlation
#' between predictors.
#'
#'Simulation framework modified from:
#'https://stats.stackexchange.com/questions/1866/how-to-simulate-a-custom-power-analysis-of-an-lm-model-using-r
#'Framework presented on stackexchange thread cites Bolker (2008) <ISBN: 9780691125220>
#'
#' @format ## `bb_type1`
#' A data frame with 76 rows and 2 columns:
#' \describe{
#'   \item{N}{Sample size}
#'   \item{Pwr}{Rejection rate}
#' }
"bb_type1"

#' Puncture: Statistical Power
#'
#' Results from puncture with a 0.5 slope of interest.
#'
#' @format ## `punc_pwr`
#' A list with two elements:
#' \describe{
#'   \item{results}{Full results from simulation}
#'   \item{summary}{Summary of results from simulation}
#' }
"punc_pwr"

#' Puncture: Type 1 Error
#'
#' Results from puncture type 1 error rate analysis
#'
#' @format ## `punc_type1`
#' A list with two elements:
#' \describe{
#'   \item{results}{Full results from simulation}
#'   \item{summary}{Summary of results from simulation}
#' }
"punc_type1"

#' Puncture: Statistical Power with Heteroskedasticity
#'
#' Results from puncture with a 0.5 slope of interest and heteroskedasticity.
#'
#' @format ## `punc_pwr_het`
#' A list with two elements:
#' \describe{
#'   \item{results}{Full results from simulation}
#'   \item{summary}{Summary of results from simulation}
#' }
"punc_pwr_het"
