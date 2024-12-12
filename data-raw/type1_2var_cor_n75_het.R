## code to prepare `type1_2var_cor_n75_het` dataset goes here
type1_2var_cor_n75_het <- puncture_lm_sim(
  n = 75,
  sim_iter = 20,
  beta_gen = function() {
    list(1,             # beta0 (intercept)
         0,             # beta1
         3)             # beta2
  },
  gen_ivs = function(n){
    mu1 <- mu2 <- 0 # means
    s1 <- 15 # sd var1
    s2 <- 10 # sd var2
    correl <- 0.6 # correlation
    Sigma <- matrix(c(s1^2, s1*s2*correl, s1*s2*correl, s2^2), 2, 2) # var-cov
    bndat <- MASS::mvrnorm(n, mu = c(mu1, mu2), Sigma = Sigma) # generate bivariate normal data
    X1 <- bndat[,1] # extract X1
    X2 <- bndat[,2] # extract X2

    #setup heteroskedastic errors
    grps = 5
    base_sd = 5
    sd_step = 10
    sequences <- split(1:n, ceiling((1:n)/(n/grps)))
    epsilon <- vector(mode = "numeric", length = n)

    # Generate increasing SDs for each group
    sds <- seq(from = base_sd, by = sd_step, length.out = grps)

    # Loop through groups to assign errors
    for(i in 1:grps) {
      group_X1 <- abs(X1[sequences[[i]]])
      epsilon[sequences[[i]]] <- stats::rnorm(n/grps, 0, group_X1*sds[i])
    }
    # epsilon <- stats::rnorm(n, 0, abs(X1))
    return(list(X1 = X1, X2 = X2, epsilon = epsilon))
  },
  b = 100,
  m = 50,
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

usethis::use_data(type1_2var_cor_n75_het, overwrite = TRUE)
