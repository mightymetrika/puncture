puncture_lm_sim <- function(n = 100,            # Sample size
                            sim_iter = 1000,      # Number of simulation iterations
                            beta_gen = function() {
                              list(1,             # beta0 (intercept)
                                   2,             # beta1
                                   3)             # beta2
                              },
                            gen_ivs = function(n) {
                              X1 <- stats::rnorm(n, 0, 1)
                              X2 <- 0.6 * X1 + sqrt(1 - 0.6^2) * stats::rnorm(n)
                              epsilon <- stats::rnorm(n)
                              list(X1 = X1, X2 = X2, epsilon = epsilon)
                              },
                            .formula = "Y ~ X1 + X2",
                            func = stats::lm,
                            combine = mean,
                            alpha = 0.05,
                            ...) {

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
    betas <- beta_gen()
    ivs <- gen_ivs(n)

    formula_vars <- all.vars(stats::as.formula(.formula))
    response_var <- formula_vars[1]
    predictor_vars <- formula_vars[-1]

    # Create data frame
    data <- as.data.frame(ivs[predictor_vars])

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
    beta1_pvalue_std <- standard_summary$coefficients[predictor_vars[1], "Pr(>|t|)"]  # New: Extract p-value

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
    # beta1_se_punct <- stats::sd(puncture_results$estimate, na.rm = TRUE)
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
