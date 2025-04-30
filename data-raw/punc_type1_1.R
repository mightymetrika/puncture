## code to prepare `punc_type1-1` dataset goes here
punc_type1_1 <- puncture_lm_sim(
  n = 20,                   # Smaller sample size (was 30)
  sim_iter = 1000,          # Keep the same number of simulation iterations
  beta_gen = function() {
    list(0,                 # beta0 (intercept)
         0,                 # beta1
         0)                 # beta2
  },
  gen_ivs = function(n){
    # Define means (all 0 for simplicity)
    mu <- rep(0, 5)         # for X1, X2, X3, X4, X5

    # Define standard deviations
    s1 <- 1                 # sd for X1 (original)
    s2 <- 3                 # sd for X2 (original)
    s3 <- 2                 # sd for X3 (auxiliary)
    s4 <- 2                 # sd for X4 (auxiliary)
    s5 <- 2                 # sd for X5 (auxiliary)

    # Define correlation matrix with mathematically valid correlations
    # Using stronger but valid correlations between variables
    R <- matrix(c(
      1.0,  0.2,  0.5,  0.4,  0.3,  # correlations with X1
      0.2,  1.0,  0.5,  0.3,  0.4,  # correlations with X2
      0.5,  0.5,  1.0,  0.2,  0.2,  # correlations with X3
      0.4,  0.3,  0.2,  1.0,  0.2,  # correlations with X4
      0.3,  0.4,  0.2,  0.2,  1.0   # correlations with X5
    ), nrow = 5, byrow = TRUE)

    # Create variance-covariance matrix
    S <- diag(c(s1, s2, s3, s4, s5))  # diagonal matrix of SDs
    Sigma <- S %*% R %*% S            # full var-covar matrix

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
      X3 = X3,
      X4 = X4,
      X5 = X5,
      epsilon = epsilon
    ))
  },
  b = 150,                  # More bootstrap iterations (was 100)
  m = 15,                   # More imputations (was 10)
  stacks = 3,               # More stacks (was 2)
  thin = c(1, 0.5, 0.25),   # Different thinning strategy
  mpat = function(mdat){    # Adjusted missingness pattern
    n <- nrow(mdat) * ncol(mdat)
    # Higher probability of missingness (lower observation probabilities)
    c1 <- stats::rbinom(n, 1, 0.75)
    c2 <- stats::rbinom(n, 1, 0.75)
    pattern_vec <- c1 * c2
    matrix(pattern_vec, nrow = nrow(mdat), ncol = ncol(mdat))
  },
  remove.collinear=FALSE,
  combine = mean
)

usethis::use_data(punc_type1_1, overwrite = TRUE)
