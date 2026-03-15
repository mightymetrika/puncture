## code to prepare `punc_recruit_type1_med.90` dataset goes here
system.time(
  punc_recruit_type1_med.90 <- puncture_lm_sim(
    n = 25,
    sim_iter = 500,
    beta_gen = function() {
      list(
        0,    # beta0 (intercept)
        0,    # beta1
        0.4   # beta2
      )
    },
    gen_ivs = function(n) {
      N_pop <- 4000L

      # Latent trait that drives the target metric and also makes recruitment harder.
      theta <- stats::rnorm(N_pop)

      # Main predictor of interest is mostly independent of theta so the example
      # focuses on compressed outcome support rather than heavy confounding.
      X1 <- stats::rnorm(N_pop)

      # Model covariate plus auxiliary variables that carry signal about theta.
      X2 <- 0.40 * theta + sqrt(1 - 0.40^2) * stats::rnorm(N_pop)
      X3 <- 0.85 * theta + sqrt(1 - 0.85^2) * stats::rnorm(N_pop)
      X4 <- 0.55 * theta + sqrt(1 - 0.55^2) * stats::rnorm(N_pop)
      X5 <- stats::rnorm(N_pop)

      # The outcome error still contains latent-trait signal, so X3 is a noisy
      # but useful predictor of Y during puncture/imputation.
      epsilon_pop <- 0.90 * theta + stats::rnorm(N_pop, 0, 0.70)

      # Higher latent values are harder to recruit, which makes larger Y values
      # less common in the observed sample on average.
      recruit_prob <- stats::plogis(0.35 - 1.10 * theta)

      idx <- sample.int(
        N_pop,
        size = n,
        replace = FALSE,
        prob = recruit_prob
      )

      list(
        X1 = X1[idx],
        X2 = X2[idx],
        X3 = X3[idx],
        X4 = X4[idx],
        X5 = X5[idx],
        epsilon = epsilon_pop[idx]
      )
    },
    b = 50,
    m = 5,
    mpat = function(mdat) {
      n <- nrow(mdat) * ncol(mdat)
      c1 <- stats::rbinom(n, 1, 0.90)
      c2 <- stats::rbinom(n, 1, 0.90)
      pattern_vec <- c1 * c2
      matrix(pattern_vec, nrow = nrow(mdat), ncol = ncol(mdat))
    },
    remove.collinear = FALSE,
    combine = median
  )
)

usethis::use_data(punc_recruit_type1_med.90, overwrite = TRUE)
