## code to prepare `type1_2var_cor_n30_sim` dataset goes here
type1_2var_cor_n30_sim <- puncture_lm_sim(
  n = 30,
  sim_iter = 10,
  beta_gen = function() {
    list(1,             # beta0 (intercept)
         0,             # beta1
         3)             # beta2
  },
  gen_ivs = function(n){
    mu1 <- mu2 <- 0 # means
    s1 <- 1 # sd var1
    s2 <- 3 # sd var2
    correl <- 0.8 # correlation
    Sigma <- matrix(c(s1^2, s1*s2*correl, s1*s2*correl, s2^2), 2, 2) # var-cov
    bndat <- MASS::mvrnorm(n, mu = c(mu1, mu2), Sigma = Sigma) # generate bivariate normal data
    X1 <- bndat[,1] # extract X1
    X2 <- bndat[,2] # extract X2
    epsilon <- stats::rnorm(n) # error term
    return(list(X1 = X1, X2 = X2, epsilon = epsilon))
  },
  b = 100,
  m = 5,
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


usethis::use_data(type1_2var_cor_n30_sim, overwrite = TRUE)
