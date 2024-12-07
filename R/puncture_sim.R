puncture_sim <- function(n = 100,            # Sample size
                         sim_iter = 1000,      # Number of simulation iterations
                         k = 5,                # Number of puncture iterations
                         b = 5,                # Number of imputations
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
                         alpha = 0.05) {

  # Initialize results storage
  results <- list(
    standard = list(
      beta1_estimates = numeric(sim_iter),
      beta1_SE = numeric(sim_iter),
      beta1_bias = numeric(sim_iter),
      beta1_CI_low = numeric(sim_iter),
      beta1_CI_high = numeric(sim_iter),
      beta1_coverage = numeric(sim_iter),
      beta1_pvalue = numeric(sim_iter),     # New: Store p-values
      beta1_rejected = numeric(sim_iter)    # New: Store rejection decisions
    ),
    puncture = list(
      beta1_estimates = numeric(sim_iter),
      beta1_SE = numeric(sim_iter),
      beta1_bias = numeric(sim_iter),
      beta1_CI_low = numeric(sim_iter),
      beta1_CI_high = numeric(sim_iter),
      beta1_coverage = numeric(sim_iter),
      beta1_pvalue = numeric(sim_iter),     # New: Store p-values
      beta1_rejected = numeric(sim_iter)
    )
  )

  # Monte Carlo simulation loop
  for (sim in 1:sim_iter) {
    # Generate data
    betas <- beta_gen()
    ivs <- gen_ivs(n)

    # Create data frame
    data <- data.frame(
      X1 = ivs$X1,
      X2 = ivs$X2
    )

    # Generate response variable
    data$Y <- betas[[1]] +
      betas[[2]] * data$X1 +
      betas[[3]] * data$X2 +
      ivs$epsilon

    # Standard analysis
    standard_model <- func(stats::as.formula(.formula), data = data)
    standard_summary <- summary(standard_model)
    beta1_hat_std <- stats::coef(standard_model)["X1"]
    beta1_se_std <- sqrt(diag(stats::vcov(standard_model)))["X1"]
    beta1_pvalue_std <- standard_summary$coefficients["X1", "Pr(>|t|)"]  # New: Extract p-value

    results$standard$beta1_estimates[sim] <- beta1_hat_std
    results$standard$beta1_SE[sim] <- beta1_se_std
    results$standard$beta1_bias[sim] <- beta1_hat_std - betas[[2]]
    results$standard$beta1_CI_low[sim] <- beta1_hat_std - stats::qnorm(1 - alpha/2) * beta1_se_std
    results$standard$beta1_CI_high[sim] <- beta1_hat_std + stats::qnorm(1 - alpha/2) * beta1_se_std
    results$standard$beta1_coverage[sim] <-
      (results$standard$beta1_CI_low[sim] <= betas[[2]]) &&
      (results$standard$beta1_CI_high[sim] >= betas[[2]])
    results$standard$beta1_pvalue[sim] <- beta1_pvalue_std               # New: Store p-value
    results$standard$beta1_rejected[sim] <- beta1_pvalue_std < alpha     # New: Store rejection decision

    # Puncture analysis
    puncture_results <- puncture(
      dat = data,
      k = k,
      b = b,
      form = .formula,
      func = func,
      term = "X1"
    )

    beta1_hat_punct <- combine(puncture_results$estimate)
    beta1_se_punct <- stats::sd(puncture_results$estimate)
    beta1_pvalue_punct <- combine(puncture_results$p.value)  # New: Combine p-values

    results$puncture$beta1_estimates[sim] <- beta1_hat_punct
    results$puncture$beta1_SE[sim] <- beta1_se_punct
    results$puncture$beta1_bias[sim] <- beta1_hat_punct - betas[[2]]
    results$puncture$beta1_CI_low[sim] <- beta1_hat_punct - stats::qnorm(1 - alpha/2) * beta1_se_punct
    results$puncture$beta1_CI_high[sim] <- beta1_hat_punct + stats::qnorm(1 - alpha/2) * beta1_se_punct
    results$puncture$beta1_coverage[sim] <-
      (results$puncture$beta1_CI_low[sim] <= betas[[2]]) &&
      (results$puncture$beta1_CI_high[sim] >= betas[[2]])
    results$puncture$beta1_pvalue[sim] <- beta1_pvalue_punct            # New: Store p-value
    results$puncture$beta1_rejected[sim] <- beta1_pvalue_punct < alpha  # New: Store rejection decision
  }

  # Compute summary statistics
  summary_stats <- data.frame(
    method = c("Standard", "Puncture"),
    avg_bias = c(
      mean(results$standard$beta1_bias),
      mean(results$puncture$beta1_bias)
    ),
    avg_se = c(
      mean(results$standard$beta1_SE),
      mean(results$puncture$beta1_SE)
    ),
    coverage = c(
      mean(results$standard$beta1_coverage),
      mean(results$puncture$beta1_coverage)
    ),
    ci_width = c(
      mean(results$standard$beta1_CI_high - results$standard$beta1_CI_low),
      mean(results$puncture$beta1_CI_high - results$puncture$beta1_CI_low)
    ),
    rejection_rate = c(                                                  # New: Add rejection rates
      mean(results$standard$beta1_rejected),
      mean(results$puncture$beta1_rejected)
    )
  )

  # Return both detailed results and summary statistics
  return(list(
    results = results,
    summary = summary_stats
  ))
}
