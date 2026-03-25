## code to prepare `punc_coarsescreener_pwr` dataset goes here

system.time(
  punc_coarsescreener_pwr <- puncture_lm_sim(
    n = 35,
    sim_iter = 500,
    beta_gen = function() {
      list(
        0,    # beta0 (intercept)
        0.35, # beta1 (target effect for power)
        0.45  # beta2
      )
    },
    gen_ivs = function(n) {
      z <- function(nn) stats::rnorm(nn)
      std <- function(x) as.numeric(scale(x))

      coarse_screener <- function(x, step = 0.5, lo = -1.5, hi = 1.5) {
        x <- pmin(pmax(x, lo), hi)
        round(x / step) * step
      }

      N_pop <- 6000L

      # Latent baseline severity / indication factor
      theta <- z(N_pop)

      # Rare complex cases that create leverage in small samples
      rare_grp <- stats::rbinom(N_pop, 1, 0.12)

      # Focal predictor is related to latent severity and has a rare-group tail
      X1 <- std(
        0.72 * theta +
          sqrt(1 - 0.72^2) * z(N_pop) +
          rare_grp * (1.40 + 0.40 * z(N_pop))
      )

      # Analyst's modeled covariate:
      # a noisy, coarse, ceiling/floor-limited baseline screener
      x2_latent <- 0.40 * theta + 0.95 * z(N_pop)
      X2 <- std(coarse_screener(x2_latent, step = 0.5, lo = -1.5, hi = 1.5))

      # Richer auxiliary intake variables: better proxies for the same severity
      X3 <- std(0.95 * theta + 0.20 * z(N_pop))
      X4 <- std(0.85 * theta + 0.30 * z(N_pop))
      X5 <- std(0.75 * theta + 0.40 * z(N_pop))
      X6 <- std(0.65 * theta + 0.50 * z(N_pop))

      # Residual severity remains in the outcome even after conditioning on coarse X2
      epsilon_pop <-
        0.60 * theta +
        stats::rnorm(N_pop, mean = 0, sd = 0.60 + 0.30 * rare_grp)

      # More severe / complex cases are somewhat harder to recruit
      recruit_prob <- stats::plogis(0.55 - 1.10 * theta - 1.15 * rare_grp)

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

      # Preserve focal signal; puncture only the nuisance/proxy block
      block_vars <- c("X2", "X3", "X4", "X5", "X6")

      # Higher-severity rows are more likely to be punctured
      block_score <- as.numeric(scale(
        0.55 * mdat$X3 +
          0.35 * mdat$X4 +
          0.20 * mdat$Y
      ))

      # Correlated blockwise puncturing
      Sigma_mask <- matrix(0.75, nrow = length(block_vars), ncol = length(block_vars))
      diag(Sigma_mask) <- 1

      latent_mask <- MASS::mvrnorm(
        n = n,
        mu = rep(0, length(block_vars)),
        Sigma = Sigma_mask
      )

      p_obs <- cbind(
        stats::plogis(0.15 - 1.05 * block_score), # X2
        stats::plogis(0.05 - 1.00 * block_score), # X3
        stats::plogis(0.20 - 0.95 * block_score), # X4
        stats::plogis(0.35 - 0.85 * block_score), # X5
        stats::plogis(0.50 - 0.75 * block_score)  # X6
      )

      p_obs <- pmax(pmin(p_obs, 0.995), 0.005)
      patt[, block_vars] <- 1L * (latent_mask < qnorm(p_obs))

      # Keep focal regressor and outcome fully observed
      patt[, "X1"] <- 1L
      patt[, "Y"] <- 1L

      patt
    },
    remove.collinear = FALSE,
    combine = mean
  )
)

usethis::use_data(punc_coarsescreener_pwr, overwrite = TRUE)
