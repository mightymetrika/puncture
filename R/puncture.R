#' Perform a Puncture Analysis
#'
#' @description
#' Implements a bootstrap-based analysis method that combines multiple imputation
#' for missing data. For each bootstrap iteration, the function creates a
#' missingness pattern ("punctures" the data), applies multiple imputation, fits a
#' specified model, and pools results across imputations.
#'
#' @param dat A data frame or matrix containing the complete dataset before introducing
#'   missing values. Must have at least 2 rows and 1 column.
#' @param b A positive integer specifying the number of bootstrap iterations to perform.
#'   Default is 5.
#' @param m A positive integer specifying the number of multiple imputations to generate
#'   for each bootstrap sample. Default is 5.
#' @param mpat A function that takes a dataset as input and returns a matrix of the same
#'   dimensions with 1's indicating observed values and 0's indicating missing values.
#'   The default function implements a simple random missingness pattern based on
#'   independent Bernoulli trials.
#' @param form A character string specifying the model formula (e.g., "y ~ x").
#' @param func A function to fit the model, such as \code{stats::lm} or \code{stats::glm}.
#'   Default is \code{stats::lm}.
#' @param term A character string specifying the predictor term of interest in the model
#'   for which statistics should be extracted.
#' @param statistics A character vector specifying which statistics to extract for the
#'   term of interest. The \code{broom::tidy()} names are used for statistics. Default
#'   is c("estimate", "std.error", "statistic", "p.value").
#' @param ... Additional arguments passed to \code{mice::mice()}.
#'
#' @return A data frame with \code{b} rows and columns corresponding to the requested
#'   statistics. Each row represents results from one bootstrap iteration after
#'   multiple imputation and pooling.
#'
#' @details
#' The function implements the following algorithm for each bootstrap iteration:
#' \enumerate{
#'   \item Generate a missingness pattern using the specified \code{mpat} function
#'   \item Apply the pattern to create a dataset with missing values
#'   \item Perform multiple imputation using \code{mice}
#'   \item Fit the specified model to each imputed dataset
#'   \item Pool results across imputations using Rubin's rules
#'   \item Extract the requested statistics for the term of interest
#' }
#'
#' @note
#' The default missingness pattern function creates missing values through two
#' independent Bernoulli trials with success probabilities 0.7 and 0.6. A value
#' is observed (1 in the pattern matrix) only if both trials are successful.
#'
#' @examples
#' puncture(cars, b = 3, m = 3, form = "speed ~ dist", term = "dist",
#'          remove.collinear=FALSE)
#'
#' @export
puncture <- function(dat,
                     b = 5,
                     m = 5,
                     mpat = function(mdat){
                       # This function creates a missingness pattern matrix
                       # based on random draws.
                       # For efficiency, consider a vectorized approach:
                       n <- nrow(mdat) * ncol(mdat)
                       # Generate vectors of random binomial draws
                       c1 <- stats::rbinom(n, 1, 0.7)
                       c2 <- stats::rbinom(n, 1, 0.6)
                       # Multiplying vectors element-wise
                       pattern_vec <- c1 * c2
                       matrix(pattern_vec, nrow = nrow(mdat), ncol = ncol(mdat))
                     },
                     form = "y ~ x",
                     func = stats::lm,
                     term = "x",
                     statistics = c("estimate", "std.error", "statistic", "p.value"),
                     ...)
{

  # Parameter checks

  ## Check dat is a data frame or matrix with rows and columns
  stopifnot("'dat' must be a data frame or matrix." = is.data.frame(dat) | is.matrix(dat))
  stopifnot("'dat' must have at least 2 rows and 1 column." = nrow(dat) >= 2 & ncol(dat) >= 1)


  ## Check k and b are positive integers
  stopifnot("'b' must be a single positive integer." =
              is.numeric(b) & length(b) == 1 & b > 0 & b == as.integer(b))
  stopifnot("'m' must be a single positive integer." =
              is.numeric(m) & length(m) == 1 & m > 0 & m == as.integer(m))

  ## Check mpat is a function
  stopifnot("'mpat' must be a function that takes 'mdat' as input and returns a missingness pattern matrix."=
              is.function(mpat))

  ## Check form is a valid formula or a character string that can be converted to a formula
  stopifnot("'form' must be a one-sided character string representing a formula."=
              is.character(form) & length(form) == 1)

  ## Check func is a function (e.g., stats::lm)
  stopifnot("'func' must be a function, for example stats::lm."=is.function(func))

  ## Check term is a single character string
  stopifnot("'term' must be a single character string representing the predictor of interest."=
              is.character(term) & length(term) == 1)

  ## Check statistics is a character vector of length >=1
  stopifnot("'statistics' must be a character vector of at least one statistic name."=
              is.character(statistics) & length(statistics) >= 1)

  # Pre-allocate result storage
  results <- data.frame(matrix(NA, nrow = b, ncol = length(statistics)))
  colnames(results) <- statistics

  # Main loop
  for (i in seq_len(b)) {
    # 1) Create missing data pattern
    temp_pattern <- mpat(dat)

    # 2) Apply pattern to data
    miss_dat <- dat
    miss_dat[temp_pattern == 0] <- NA

    # 3) Perform multiple imputation
    midat <- mice::mice(miss_dat, m = m, printFlag = FALSE, ...)

    # 4) Fit model to each imputed dataset, then pool results
    mimod <- with(midat, func((stats::as.formula(form))))
    pooled <- mice::pool(mimod)
    pooled_tidy <- broom::tidy(pooled)

    # 5) Extract requested statistics for specified term
    term_row <- pooled_tidy[pooled_tidy$term == term, statistics, drop = FALSE]
    results[i, ] <- as.vector(term_row)
  }

  # Return only the bootstrap replicates as requested
  return(results)
}
