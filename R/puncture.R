puncture <- function(dat,
                     k = 5,
                     b = 5,
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
  stopifnot("'k' must be a single positive integer." =
              is.numeric(k) & length(k) == 1 & k > 0 & k == as.integer(k))
  stopifnot("'b' must be a single positive integer." =
              is.numeric(b) & length(b) == 1 & b > 0 & b == as.integer(b))

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
  results <- data.frame(matrix(NA, nrow = k, ncol = length(statistics)))
  colnames(results) <- statistics

  # Main loop
  for (i in seq_len(k)) {
    # 1) Create missing data pattern
    temp_pattern <- mpat(dat)

    # 2) Apply pattern to data
    miss_dat <- dat
    miss_dat[temp_pattern == 0] <- NA

    # 3) Perform multiple imputation
    midat <- mice::mice(miss_dat, m = b, printFlag = FALSE, ...)

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
