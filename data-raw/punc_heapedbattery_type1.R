## code to prepare `punc_heapedbattery_type1` dataset goes here

system.time(
  punc_heapedbattery_type1 <- puncture_lm_sim(
    n = 35,
    sim_iter = 500,
    beta_gen = function() {
      list(
        0,    # beta0 (intercept)
        0,    # beta1 (null effect for Type I)
        0.50  # beta2
      )
    },
    gen_ivs = function(n) {
      z <- function(nn) stats::rnorm(nn)
      std <- function(x) as.numeric(scale(x))

      # Short self-report screener:
      # upper-tail compression + heaping on round values + mild top-coding.
      heap_tail <- function(x) {
        x_obs <- x

        # Compress the upper tail (e.g., ceiling / under-reporting among severe cases)
        hi <- x > 1.10
        x_obs[hi] <- 1.10 + 0.30 * (x[hi] - 1.10)

        # More heaping in the upper tail
        heap_prob <- stats::plogis(-0.40 + 1.20 * pmax(x_obs, 0))
        heaped <- round(x_obs * 2) / 2

        use_heap <- stats::rbinom(length(x_obs), 1, heap_prob)
        x_mix <- ifelse(
          use_heap == 1,
          heaped,
          x_obs + stats::rnorm(length(x_obs), mean = 0, sd = 0.08)
        )

        # Mild top-code
        pmin(pmax(x_mix, -2.50), 1.50)
      }

      N_pop <- 7000L

      # Latent baseline severity / indication factor
      theta <- z(N_pop)

      # Rare complex subgroup that is harder to recruit and carries extra leverage
      rare_grp <- stats::rbinom(N_pop, 1, 0.10)

      # Exposure / treatment dose / uptake:
      # correlated with latent severity (confounding by indication)
      X1 <- std(
        0.68 * theta +
          sqrt(1 - 0.68^2) * z(N_pop) +
          rare_grp * (1.25 + 0.35 * z(N_pop))
      )

      # Analyst's modeled covariate:
      # a short self-report screener measured worst in the upper tail
      x2_latent <- 0.90 * theta + 0.20 * z(N_pop) + 0.15 * rare_grp
      X2 <- std(heap_tail(x2_latent))

      # Better auxiliary proxies for the same latent severity process
      X3 <- std(0.95 * theta + 0.18 * z(N_pop)) # clinician rating
      X4 <- std(0.88 * theta + 0.25 * z(N_pop)) # caregiver / collateral rating
      X5 <- std(0.80 * theta + 0.35 * z(N_pop)) # passive sensor / behavioral trace
      X6 <- std(0.72 * theta + 0.45 * z(N_pop)) # biomarker / functional measure

      # Residual severity remains after conditioning on distorted X2
      epsilon_pop <-
        0.55 * theta +
        stats::rnorm(N_pop, mean = 0, sd = 0.60 + 0.30 * rare_grp)

      # Selective recruitment compresses the upper tail of severity in the observed sample
      recruit_prob <- stats::plogis(
        0.55 - 1.20 * pmax(theta, 0) - 1.10 * rare_grp
      )

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

      # Preserve focal signal; puncture only the confounder/proxy block
      block_vars <- c("X2", "X3", "X4", "X5", "X6")

      # Extra puncturing for rows likely to be in the distorted upper tail
      top_code_flag <- as.numeric(
        mdat$X2 > stats::quantile(mdat$X2, probs = 0.80, na.rm = TRUE)
      )

      block_score <- as.numeric(scale(
        0.45 * mdat$X3 +
          0.30 * mdat$X4 +
          0.20 * mdat$Y +
          0.35 * top_code_flag
      ))

      # Correlated blockwise puncturing
      Sigma_mask <- matrix(
        0.70,
        nrow = length(block_vars),
        ncol = length(block_vars)
      )
      diag(Sigma_mask) <- 1

      latent_mask <- MASS::mvrnorm(
        n = n,
        mu = rep(0, length(block_vars)),
        Sigma = Sigma_mask
      )

      p_obs <- cbind(
        stats::plogis(0.10 - 1.30 * block_score), # X2: puncture hardest
        stats::plogis(0.15 - 1.10 * block_score), # X3
        stats::plogis(0.25 - 1.00 * block_score), # X4
        stats::plogis(0.40 - 0.90 * block_score), # X5
        stats::plogis(0.50 - 0.80 * block_score)  # X6
      )

      p_obs <- pmax(pmin(p_obs, 0.995), 0.005)
      patt[, block_vars] <- 1L * (latent_mask < stats::qnorm(p_obs))

      # Keep target outcome and focal regressor fully observed
      patt[, "X1"] <- 1L
      patt[, "Y"] <- 1L

      patt
    },
    remove.collinear = FALSE,
    combine = mean
  )
)

usethis::use_data(punc_heapedbattery_type1, overwrite = TRUE)
