# data-raw/punc_latentproxy_pwr.R

## code to prepare `punc_latentproxy_pwr` dataset goes here

system.time(
  punc_latentproxy_pwr <- puncture_lm_sim(
    n = 25,
    sim_iter = 500,
    beta_gen = function() {
      list(
        0,    # beta0 (intercept)
        0.45, # beta1 (target effect for power)
        0.45  # beta2
      )
    },
    gen_ivs = function(n) {
      z <- function(nn) stats::rnorm(nn)
      std <- function(x) as.numeric(scale(x))

      N_pop <- 6000L

      # Latent severity / nuisance trait
      theta <- z(N_pop)

      # Rare complex cases that create extra leverage in small samples
      rare_grp <- stats::rbinom(N_pop, 1, 0.10)

      # Focal predictor is correlated with latent severity and has a rare-group tail.
      X1 <- std(
        0.70 * theta +
          sqrt(1 - 0.70^2) * z(N_pop) +
          rare_grp * (1.60 + 0.45 * z(N_pop))
      )

      # Analyst's modeled covariate: intentionally weak / noisy proxy for theta.
      X2 <- std(0.20 * theta + 0.98 * z(N_pop))

      # Auxiliary variables: much stronger proxies for the same latent severity.
      X3 <- std(0.95 * theta + 0.20 * z(N_pop))
      X4 <- std(0.85 * theta + 0.30 * z(N_pop))
      X5 <- std(0.75 * theta + 0.40 * z(N_pop))
      X6 <- std(0.65 * theta + 0.50 * z(N_pop))

      # Residual still carries latent severity signal even after conditioning on X2.
      epsilon_pop <-
        0.70 * theta +
        stats::rnorm(N_pop, mean = 0, sd = 0.65 + 0.30 * rare_grp)

      # More severe / rare-group units are somewhat harder to recruit,
      # compressing support in the realized small sample.
      recruit_prob <- stats::plogis(0.40 - 1.00 * theta - 0.80 * rare_grp)

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

      # Preserve the focal signal and puncture only the nuisance/proxy block.
      block_vars <- c("X2", "X3", "X4", "X5", "X6")

      # Higher-severity rows are more likely to be punctured.
      block_score <- as.numeric(scale(
        0.60 * mdat$X3 +
          0.40 * mdat$X4 +
          0.20 * mdat$Y
      ))

      # Correlated blockwise puncturing across the nuisance/proxy variables.
      Sigma_mask <- matrix(0.75, nrow = length(block_vars), ncol = length(block_vars))
      diag(Sigma_mask) <- 1

      latent_mask <- MASS::mvrnorm(
        n = n,
        mu = rep(0, length(block_vars)),
        Sigma = Sigma_mask
      )

      p_obs <- cbind(
        stats::plogis(0.00 - 1.20 * block_score),  # X2
        stats::plogis(0.10 - 1.10 * block_score),  # X3
        stats::plogis(0.25 - 1.00 * block_score),  # X4
        stats::plogis(0.40 - 0.90 * block_score),  # X5
        stats::plogis(0.55 - 0.80 * block_score)   # X6
      )

      p_obs <- pmax(pmin(p_obs, 0.995), 0.005)
      patt[, block_vars] <- 1L * (latent_mask < qnorm(p_obs))

      # Keep the target outcome and focal regressor fully observed.
      patt[, "X1"] <- 1L
      patt[, "Y"] <- 1L

      patt
    },
    remove.collinear = FALSE,
    combine = mean
  )
)

usethis::use_data(punc_latentproxy_pwr, overwrite = TRUE)
