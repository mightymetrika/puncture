## code to prepare `type1_2var_uncor_n30_het` dataset goes here
type1_2var_uncor_n30_het <- puncture_lm_sim(
  n = 30,
  sim_iter = 10,
  beta_gen = function() {
    list(1,             # beta0 (intercept)
         0,             # beta1
         3)             # beta2
  },
  gen_ivs = function(n){
    X1 <- stats::rnorm(n, mean = 0, sd = 5)
    X2 <- stats::rnorm(n, 0, 3)
    epsilon <- stats::rnorm(n, 0, abs(X1)) # error term
    return(list(X1 = X1, X2 = X2, epsilon = epsilon))
  },
  b = 100,
  m = 40,
  mpat = function(mdat){ # define missing data pattern
    n <- nrow(mdat) * ncol(mdat)
    c1 <- stats::rbinom(n, 1, 0.8)
    c2 <- stats::rbinom(n, 1, 0.8)
    pattern_vec <- c1 * c2
    matrix(pattern_vec, nrow = nrow(mdat), ncol = ncol(mdat))
  },
  remove.collinear=FALSE,
  combine = mean
)

usethis::use_data(type1_2var_uncor_n30_het, overwrite = TRUE)
