## code to prepare `punc_proxyseverity_type1` dataset goes here

system.time(
  punc_proxyseverity_type1 <- puncture_lm_sim(
    n = 25,
    sim_iter = 500,
    beta_gen = function() {
      list(
        0,    # beta0 (intercept)
        0,    # beta1
        0.35  # beta2
      )
    },
    gen_ivs = function(n) {
      z <- function(nn) stats::rnorm(nn)
      std <- function(x) as.numeric(scale(x))

      N_pop <- 5000L

      # Latent baseline severity that is only partially measured by the modeled covariate.
      theta <- z(N_pop)

      # Rare complex cases create extra leverage in small samples and are somewhat
      # harder to recruit into observational pilot studies.
      rare_grp <- stats::rbinom(N_pop, 1, 0.10)

      # Focal predictor is related to underlying severity and has a small rare-group tail.
      X1 <- std(
        0.65 * theta +
          sqrt(1 - 0.65^2) * z(N_pop) +
          rare_grp * (1.80 + 0.50 * z(N_pop))
      )

      # X2 is the analyst's short baseline screener: useful, but noisy.
      X2 <- std(0.35 * theta + 0.90 * z(N_pop))

      # Auxiliary variables are stronger proxies for the same latent severity process.
      X3 <- std(0.92 * theta + 0.25 * z(N_pop))
      X4 <- std(0.82 * theta + 0.35 * z(N_pop))
      X5 <- std(0.72 * theta + 0.45 * z(N_pop))
      X6 <- std(0.60 * theta + 0.55 * z(N_pop))

      # Residual severity signal remains in the outcome even after controlling for X2,
      # which makes the single noisy screener an imperfect adjustment set.
      epsilon_pop <-
        0.55 * theta +
        stats::rnorm(N_pop, mean = 0, sd = 0.70 + 0.35 * rare_grp)

      # Higher-severity and rare-group units are less likely to be recruited,
      # compressing support in the realized sample.
      recruit_prob <- stats::plogis(0.45 - 1.10 * theta - 1.00 * rare_grp)

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
        X6 = X6[idx],
        epsilon = epsilon_pop[idx]
      )
    },
    b = 50,
    m = 5,
    mpat = function(mdat) {
      n <- nrow(mdat)
      p <- ncol(mdat)

      patt <- matrix(1L, nrow = n, ncol = p)
      colnames(patt) <- colnames(mdat)

      # Preserve the target signal and puncture only the nuisance/proxy block.
      block_vars <- c("X2", "X3", "X4", "X5", "X6")

      # Severity-related rows are more likely to be punctured, creating a blockwise
      # missingness structure that multiple auxiliaries can help repair.
      block_score <- as.numeric(scale(
        0.50 * mdat$X3 +
          0.35 * mdat$X4 +
          0.20 * mdat$Y
      ))

      Sigma_mask <- matrix(0.70, nrow = length(block_vars), ncol = length(block_vars))
      diag(Sigma_mask) <- 1
      latent_mask <- MASS::mvrnorm(
        n = n,
        mu = rep(0, length(block_vars)),
        Sigma = Sigma_mask
      )

      p_obs <- cbind(
        stats::plogis(0.15 - 1.10 * block_score),
        stats::plogis(0.05 - 1.05 * block_score),
        stats::plogis(0.20 - 0.95 * block_score),
        stats::plogis(0.35 - 0.90 * block_score),
        stats::plogis(0.45 - 0.80 * block_score)
      )

      p_obs <- pmax(pmin(p_obs, 0.995), 0.005)
      patt[, block_vars] <- 1L * (latent_mask < qnorm(p_obs))

      patt[, "X1"] <- 1L
      patt[, "Y"] <- 1L

      patt
    },
    remove.collinear = FALSE,
    combine = mean
  )
)

usethis::use_data(punc_proxyseverity_type1, overwrite = TRUE)

