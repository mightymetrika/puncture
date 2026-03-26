## code to prepare `punc_tailcalib_type1` dataset goes here

std <- function(x) as.numeric(scale(x))

coarse_tail_screener <- function(x, step = 0.75, lo = -2.25, hi = 2.25) {
  x <- pmin(pmax(x, lo), hi)
  round(x / step) * step
}

make_mice_spec <- function() {
  vars <- c("X1", "X2", "X3", "X4", "X5", "X6", "Y")

  method <- rep("", length(vars))
  names(method) <- vars
  method["X2"] <- "pmm"

  predictorMatrix <- matrix(
    0,
    nrow = length(vars),
    ncol = length(vars),
    dimnames = list(vars, vars)
  )
  predictorMatrix["X2", c("X1", "X3", "X4", "X5", "X6", "Y")] <- 1
  diag(predictorMatrix) <- 0

  list(
    method = method,
    predictorMatrix = predictorMatrix
  )
}

mice_spec <- make_mice_spec()

system.time(
  punc_tailcalib_type1 <- puncture_lm_sim(
    n = 30,
    sim_iter = 500,
    beta_gen = function() {
      list(
        0.00, # beta0
        0.00, # beta1
        0.60  # beta2
      )
    },
    gen_ivs = function(n) {
      z <- function(nn) stats::rnorm(nn)

      N_pop <- 10000L

      # Rare lower-tail and upper-tail subgroups.
      tail_sign <- sample(
        c(-1, 0, 1),
        size = N_pop,
        replace = TRUE,
        prob = c(0.10, 0.80, 0.10)
      )

      theta <- z(N_pop) + tail_sign * 2.35
      tail_flag <- as.integer(tail_sign != 0)

      # Focal predictor correlated with latent risk.
      X1 <- std(
        0.60 * theta +
          sqrt(1 - 0.60^2) * z(N_pop) +
          0.35 * tail_flag * z(N_pop)
      )

      # Analyst's actual adjustment covariate:
      # short, coarse, saturated screener that compresses the tails.
      x2_latent <- 0.35 * theta +
        0.30 * sign(theta) * pmax(abs(theta) - 1.00, 0) +
        0.95 * z(N_pop)

      X2 <- std(
        coarse_tail_screener(
          x2_latent,
          step = 0.75,
          lo = -2.25,
          hi = 2.25
        )
      )

      # Rich intake battery / auxiliary variables that preserve tail information.
      X3 <- std(theta + 0.15 * z(N_pop))
      X4 <- std(0.90 * theta + 0.25 * z(N_pop))
      X5 <- std(sign(theta) * abs(theta)^0.80 + 0.30 * z(N_pop))
      X6 <- std(0.70 * theta + 0.45 * z(N_pop))

      # Outcome still depends on latent tail structure beyond what coarse X2 captures.
      epsilon_pop <-
        0.45 * theta +
        0.35 * (abs(theta) - mean(abs(theta))) +
        stats::rnorm(
          N_pop,
          mean = 0,
          sd = 0.55 + 0.20 * tail_flag
        )

      # Both tails are harder to recruit, so small samples underrepresent them.
      recruit_prob <- stats::plogis(
        1.00 - 0.95 * abs(theta) - 0.55 * tail_flag
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
    b = 80,
    m = 10,
    mpat = function(mdat) {
      patt <- matrix(1L, nrow = nrow(mdat), ncol = ncol(mdat))
      colnames(patt) <- colnames(mdat)

      # Tail-like rows are the ones where the coarse screener is least trustworthy.
      tail_score <- as.numeric(scale(
        0.45 * abs(mdat$X3) +
          0.30 * abs(mdat$X4) +
          0.15 * abs(mdat$X5) +
          0.10 * abs(mdat$Y)
      ))

      # Puncture only X2, not the auxiliaries.
      p_obs_x2 <- stats::plogis(0.65 - 1.60 * tail_score)
      p_obs_x2 <- pmax(pmin(p_obs_x2, 0.995), 0.005)

      patt[, "X2"] <- stats::rbinom(nrow(mdat), 1, p_obs_x2)
      patt[, setdiff(colnames(patt), "X2")] <- 1L

      patt
    },
    remove.collinear = FALSE,
    combine = median,
    method = mice_spec$method,
    predictorMatrix = mice_spec$predictorMatrix,
    maxit = 10
  )
)

usethis::use_data(punc_tailcalib_type1, overwrite = TRUE)
