#' Simulation Framework for Comparing Standard and Puncture Methods
#'
#' @description
#' Conducts Monte Carlo simulations to compare the performance of standard linear
#' regression against the puncture method. Generates data under specified conditions,
#' applies both methods, and computes various performance metrics including bias,
#' coverage, and rejection rates.
#'
#' @param n Positive integer specifying the sample size for each simulation iteration.
#'   Default is 100.
#' @param sim_iter Positive integer specifying the number of Monte Carlo iterations.
#'   Default is 1000.
#' @param beta_gen Function that returns a list of true population parameters:
#'   first element is the intercept (beta0), followed by slope coefficients
#'   (beta1, beta2, etc.).
#' @param gen_ivs Function that takes sample size n as input and returns a list
#'   containing the predictor variables and error term (epsilon).
#' @param .formula Character string specifying the model formula. Default is
#'   "Y ~ X1 + X2".
#' @param func Function used to fit the model (e.g., stats::lm). Default is
#'   stats::lm.
#' @param combine Function used to combine results from multiple bootstrap
#'   iterations in the puncture method. Default is mean.
#' @param alpha Numeric value between 0 and 1 specifying the significance level
#'   for hypothesis tests and confidence intervals. Default is 0.05.
#' @param ... Additional arguments passed to the puncture function.
#'
#' @return A list containing two elements:
#' \describe{
#'   \item{results}{A list with detailed results for both methods, including:
#'     \itemize{
#'       \item beta1_estimates: Point estimates
#'       \item beta1_SE: Standard errors
#'       \item beta1_bias: Bias from true parameter
#'       \item beta1_CI_low: Lower confidence interval bounds
#'       \item beta1_CI_high: Upper confidence interval bounds
#'       \item beta1_coverage: Indicator for CI coverage of true parameter
#'       \item beta1_pvalue: P-values for hypothesis tests
#'       \item beta1_rejected: Indicator for null hypothesis rejection
#'       \item converged: Convergence indicators
#'     }
#'   }
#'   \item{summary}{A data frame comparing both methods on:
#'     \itemize{
#'       \item Average bias
#'       \item Average standard error
#'       \item Coverage probability
#'       \item Average confidence interval width
#'       \item Rejection rate
#'       \item Convergence rate
#'     }
#'   }
#' }
#'
#' @details
#' The function implements the following simulation procedure:
#' \enumerate{
#'   \item Generates true parameters using beta_gen()
#'   \item Generates predictor variables and error term using gen_ivs()
#'   \item Constructs response variable using linear model
#'   \item Fits standard linear model
#'   \item Applies puncture method
#'   \item Computes and stores performance metrics
#'   \item Repeats for specified number of iterations
#'   \item Summarizes results across all iterations
#' }
#'
#' @examples
#' sim_results <- puncture_lm_sim(n = 200, sim_iter = 5,
#' beta_gen = function() {
#'   list(1,             # beta0 (intercept)
#'        1,             # beta1
#'        3)             # beta2
#' },
#' gen_ivs = function(n) {
#'   X1 <- stats::rnorm(n, 0, 20)
#'   X2 <- 0.6 * X1 + sqrt(1 - 0.6^2) * stats::rnorm(n)
#'   epsilon <- stats::rnorm(n, 0, abs(X1))
#'   list(X1 = X1, X2 = X2, epsilon = epsilon)
#' },
#' m = 3, b = 3) |> suppressWarnings()
#'
#' @export
puncture_lm_sim <- function(n = 100,            # Sample size
                            sim_iter = 1000,      # Number of simulation iterations
                            beta_gen,
                            gen_ivs,
                            .formula = "Y ~ X1 + X2",
                            func = stats::lm,
                            combine = mean,
                            alpha = 0.05,
                            ...) {

  # Parameter checks

  ## Check n is a positive integer
  stopifnot("'n' must be a single positive integer." =
              is.numeric(n) && length(n) == 1 && n > 0 && n == as.integer(n))

  ## Check sim_iter is a positive integer
  stopifnot("'sim_iter' must be a single positive integer." =
              is.numeric(sim_iter) && length(sim_iter) == 1 && sim_iter > 0 &&
              sim_iter == as.integer(sim_iter))

  ## Check beta_gen is a function
  stopifnot("'beta_gen' must be a function." = is.function(beta_gen))

  ## Check gen_ivs is a function
  stopifnot("'gen_ivs' must be a function." = is.function(gen_ivs))

  ## Check func is a function
  stopifnot("'func' must be a function (e.g. stats::lm)." = is.function(func))

  ## Check combine is a function
  stopifnot("'combine' must be a function." = is.function(combine))

  ## Check alpha is numeric, length 1, and in (0,1)
  stopifnot("'alpha' must be a single numeric value between 0 and 1." =
              is.numeric(alpha) && length(alpha) == 1 && alpha > 0 && alpha < 1)

  # Initialize results storage
  results <- list(
    standard = list(
      beta1_estimates = numeric(sim_iter),
      beta1_SE = numeric(sim_iter),
      beta1_bias = numeric(sim_iter),
      beta1_CI_low = numeric(sim_iter),
      beta1_CI_high = numeric(sim_iter),
      beta1_coverage = numeric(sim_iter),
      beta1_pvalue = numeric(sim_iter),
      beta1_rejected = numeric(sim_iter),
      converged = numeric(sim_iter)
    ),
    puncture = list(
      beta1_estimates = numeric(sim_iter),
      beta1_SE = numeric(sim_iter),
      beta1_bias = numeric(sim_iter),
      beta1_CI_low = numeric(sim_iter),
      beta1_CI_high = numeric(sim_iter),
      beta1_coverage = numeric(sim_iter),
      beta1_pvalue = numeric(sim_iter),
      beta1_rejected = numeric(sim_iter),
      converged = numeric(sim_iter)
    )
  )

  # Monte Carlo simulation loop
  for (sim in 1:sim_iter) {
    # Generate data
    betas <- beta_gen() # get true population parameters
    if (!is.list(betas) || length(betas) < 2) {
      stop("'beta_gen' must return a list of at least two elements: intercept and at least one slope.")
    }

    ivs <- gen_ivs(n) # generate independent variables
    ## Check that gen_ivs returns a list and that epsilon is included
    if (!is.list(ivs) || !"epsilon" %in% names(ivs)) {
      stop("'gen_ivs' must return a list containing 'epsilon'.")
    }

    formula_vars <- all.vars(stats::as.formula(.formula))
    response_var <- formula_vars[1]
    predictor_vars <- formula_vars[-1]

    # Create data frame
    data <- as.data.frame(ivs[predictor_vars])

    ## Check that the required predictor variables exist
    if (!all(predictor_vars %in% names(data))) {
      stop("The predictor variables in the formula are not all available in the data generated by 'gen_ivs'.")
    }

    # Generate response variable
    data[[response_var]] <- 0
    for (i in seq_along(betas)) {
      if (i == 1) {
        data[[response_var]] <- data[[response_var]] + betas[[i]]  # Intercept
      } else {
        data[[response_var]] <- data[[response_var]] + betas[[i]] * data[[predictor_vars[i-1]]]
      }
    }
    data[[response_var]] <- data[[response_var]] + ivs$epsilon  # Add error term

    # Standard analysis
    standard_model <- stats::lm(stats::as.formula(.formula), data = data)
    standard_summary <- summary(standard_model)
    beta1_hat_std <- stats::coef(standard_model)[predictor_vars[1]]
    beta1_se_std <- sqrt(diag(stats::vcov(standard_model)))[[predictor_vars[1]]]
    beta1_pvalue_std <- standard_summary$coefficients[predictor_vars[1], "Pr(>|t|)"]

    results$standard$beta1_estimates[sim] <- beta1_hat_std
    results$standard$beta1_SE[sim] <- beta1_se_std
    results$standard$beta1_bias[sim] <- beta1_hat_std - betas[[2]]
    results$standard$beta1_CI_low[sim] <- beta1_hat_std - stats::qnorm(1 - alpha/2) * beta1_se_std
    results$standard$beta1_CI_high[sim] <- beta1_hat_std + stats::qnorm(1 - alpha/2) * beta1_se_std
    results$standard$beta1_coverage[sim] <-
      (results$standard$beta1_CI_low[sim] <= betas[[2]]) &&
      (results$standard$beta1_CI_high[sim] >= betas[[2]])
    results$standard$beta1_pvalue[sim] <- beta1_pvalue_std
    results$standard$beta1_rejected[sim] <- beta1_pvalue_std < alpha
    results$standard$converged[sim] <- ifelse(!is.na(beta1_hat_std), 1, 0)

    # Puncture analysis
    puncture_results <- puncture(
      dat = data,
      form = .formula,
      func = stats::lm,
      term = predictor_vars[1],
      ...
    )

    beta1_hat_punct <- combine(puncture_results$estimate, na.rm = TRUE)
    beta1_se_punct <- combine(puncture_results$std.error, na.rm = TRUE)
    beta1_pvalue_punct <- combine(puncture_results$p.value, na.rm = TRUE)

    results$puncture$beta1_estimates[sim] <- beta1_hat_punct
    results$puncture$beta1_SE[sim] <- beta1_se_punct
    results$puncture$beta1_bias[sim] <- beta1_hat_punct - betas[[2]]
    results$puncture$beta1_CI_low[sim] <- beta1_hat_punct - stats::qnorm(1 - alpha/2) * beta1_se_punct
    results$puncture$beta1_CI_high[sim] <- beta1_hat_punct + stats::qnorm(1 - alpha/2) * beta1_se_punct
    results$puncture$beta1_coverage[sim] <-
      (results$puncture$beta1_CI_low[sim] <= betas[[2]]) &&
      (results$puncture$beta1_CI_high[sim] >= betas[[2]])
    results$puncture$beta1_pvalue[sim] <- beta1_pvalue_punct
    results$puncture$beta1_rejected[sim] <- beta1_pvalue_punct < alpha
    results$puncture$converged[sim] <- mean(!is.na(puncture_results$estimate), na.rm = TRUE)
  }

  # Compute summary statistics
  summary_stats <- data.frame(
    method = c("Standard", "Puncture"),
    avg_bias = c(
      mean(results$standard$beta1_bias, na.rm = TRUE),
      mean(results$puncture$beta1_bias, na.rm = TRUE)
    ),
    avg_se = c(
      mean(results$standard$beta1_SE, na.rm = TRUE),
      mean(results$puncture$beta1_SE, na.rm = TRUE)
    ),
    coverage = c(
      mean(results$standard$beta1_coverage, na.rm = TRUE),
      mean(results$puncture$beta1_coverage, na.rm = TRUE)
    ),
    ci_width = c(
      mean(results$standard$beta1_CI_high - results$standard$beta1_CI_low, na.rm = TRUE),
      mean(results$puncture$beta1_CI_high - results$puncture$beta1_CI_low, na.rm = TRUE)
    ),
    rejection_rate = c(
      mean(results$standard$beta1_rejected, na.rm = TRUE),
      mean(results$puncture$beta1_rejected, na.rm = TRUE)
    ),
    converged = c(
      mean(results$standard$converged, na.rm = TRUE),
      mean(results$puncture$converged, na.rm = TRUE)
    )
  )

  # Return both detailed results and summary statistics
  return(list(
    results = results,
    summary = summary_stats
  ))
}
