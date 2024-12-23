## code to prepare `punc_pwr` dataset goes here
system.time(
punc_pwr <- puncture_lm_sim(
  n = 30,
  sim_iter = 1000,
  beta_gen = function() {
    list(0,             # beta0 (intercept)
         0.5,             # beta1
         0)             # beta2
  },
  gen_ivs = function(n){
    # Define means (all 0 for simplicity)
    mu <- rep(0, 5)  # for X1, X2, X3, X4, X5

    # Define standard deviations
    s1 <- 1  # sd for X1 (original)
    s2 <- 3  # sd for X2 (original)
    s3 <- 2  # sd for X3 (auxiliary)
    s4 <- 2  # sd for X4 (auxiliary)
    s5 <- 2  # sd for X5 (auxiliary)

    # Define correlation matrix
    # Using moderate correlations between auxiliary vars and main predictors
    R <- matrix(c(
      1.0,  0.1,  0.3,  0.4,  0.2,  # correlations with X1
      0.1,  1.0,  0.4,  0.2,  0.3,  # correlations with X2
      0.3,  0.4,  1.0,  0.1,  0.1,  # correlations with X3
      0.4,  0.2,  0.1,  1.0,  0.1,  # correlations with X4
      0.2,  0.3,  0.1,  0.1,  1.0   # correlations with X5
    ), nrow = 5, byrow = TRUE)

    # Create variance-covariance matrix
    S <- diag(c(s1, s2, s3, s4, s5))  # diagonal matrix of SDs
    Sigma <- S %*% R %*% S  # full var-covar matrix

    # Generate multivariate normal data
    bndat <- MASS::mvrnorm(n, mu = mu, Sigma = Sigma)

    # Extract variables
    X1 <- bndat[,1]  # main predictor 1
    X2 <- bndat[,2]  # main predictor 2
    X3 <- bndat[,3]  # auxiliary variable
    X4 <- bndat[,4]  # auxiliary variable
    X5 <- bndat[,5]  # auxiliary variable

    # Generate error term
    epsilon <- stats::rnorm(n)

    # Return all variables
    return(list(
      X1 = X1,
      X2 = X2,
      X3 = X3,  # including auxiliary variables
      X4 = X4,
      X5 = X5,
      epsilon = epsilon
    ))
  },
  b = 1000,
  m = 40,
  stacks = 2,
  thin = c(1, 1/3),
  mpat = function(mdat){ # define missing data pattern
    n <- nrow(mdat) * ncol(mdat)
    c1 <- stats::rbinom(n, 1, 0.77)
    c2 <- stats::rbinom(n, 1, 0.78)
    pattern_vec <- c1 * c2
    matrix(pattern_vec, nrow = nrow(mdat), ncol = ncol(mdat))
  },
  remove.collinear=FALSE,
  combine = mean
)
)
usethis::use_data(punc_pwr, overwrite = TRUE)
