## code to prepare `punc_anchorblock_exoerr_pwr` dataset goes here
system.time(
  punc_anchorblock_exoerr_pwr <- puncture_lm_sim(
    n = 25,
    sim_iter = 500,
    beta_gen = function() {
      list(
        0,    # beta0 (intercept)
        0.50, # beta1
        0.40  # beta2
      )
    },
    gen_ivs = function(n) {
      z <- function(nn) stats::rnorm(nn)
      std <- function(x) as.numeric(scale(x))

      N_pop <- 6000L

      # Latent nuisance factor that the auxiliary variables measure well.
      theta <- z(N_pop)

      # Rare high-leverage subgroup. These units carry unusually extreme X1 values,
      # but are harder to recruit, so the realized sample often contains only a few.
      rare_grp <- stats::rbinom(N_pop, 1, 0.12)

      # Focal regressor with a rare upper-tail subgroup.
      x1_base <- 0.70 * theta + sqrt(1 - 0.70^2) * z(N_pop)
      X1 <- std(x1_base + rare_grp * (2.75 + 0.90 * z(N_pop)))

      # Strong near-collinearity in the modeled regressors.
      X2 <- std(0.92 * X1 + 0.25 * theta + 0.20 * z(N_pop))

      # Anchor auxiliaries: strong measurements of theta and X2.
      X3 <- std(0.90 * theta + 0.25 * X2 + 0.25 * z(N_pop))
      X4 <- std(0.80 * theta + 0.35 * X2 + 0.35 * z(N_pop))
      X5 <- std(0.65 * theta + 0.50 * X2 + 0.45 * z(N_pop))
      X6 <- std(0.55 * theta + 0.60 * X2 + 0.55 * z(N_pop))

      # Outcome noise still contains latent-factor structure, so auxiliaries are useful
      # during imputation when the nuisance block is punctured.
      # Still difficult structure like punc_anchorblock_pwr, but error term is exogenous
      epsilon_pop <- rnorm(N_pop, mean = 0, sd = 0.60 + 0.50 * rare_grp)

      # Recruitment under-samples the rare leverage cases and mildly under-samples
      # higher-theta units.
      recruit_prob <- stats::plogis(0.25 - 1.40 * rare_grp - 0.35 * theta)

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

      # Best-shot mask for puncture in the current framework:
      # keep Y and X1 fully observed, then puncture only the nuisance block.
      block_vars <- c("X2", "X3", "X4", "X5", "X6")

      # Shared row-level signal for informative missingness.
      block_score <- as.numeric(scale(
        0.55 * mdat$X3 +
          0.40 * mdat$X4 +
          0.25 * mdat$X5 +
          0.15 * mdat$Y
      ))

      # Correlated latent Gaussian mask to induce blockwise missingness.
      Sigma_mask <- matrix(0.65, nrow = length(block_vars), ncol = length(block_vars))
      diag(Sigma_mask) <- 1
      latent_mask <- MASS::mvrnorm(
        n = n,
        mu = rep(0, length(block_vars)),
        Sigma = Sigma_mask
      )

      p_obs <- cbind(
        stats::plogis(0.40 - 1.10 * block_score),
        stats::plogis(0.10 - 1.15 * block_score),
        stats::plogis(0.20 - 1.00 * block_score),
        stats::plogis(0.35 - 0.95 * block_score),
        stats::plogis(0.45 - 0.85 * block_score)
      )

      p_obs <- pmax(pmin(p_obs, 0.995), 0.005)
      patt[, block_vars] <- 1L * (latent_mask < qnorm(p_obs))

      # Explicitly preserve the target signal columns.
      patt[, "X1"] <- 1L
      patt[, "Y"] <- 1L

      patt
    },
    remove.collinear = FALSE,
    combine = mean
  )
)

usethis::use_data(punc_anchorblock_exoerr_pwr, overwrite = TRUE)
