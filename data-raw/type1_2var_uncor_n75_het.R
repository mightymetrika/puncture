## code to prepare `type1_2var_uncor_n75_het` dataset goes here
type1_2var_uncor_n75_het <- puncture_lm_sim(
  n = 75,
  sim_iter = 10,
  beta_gen = function() {
    list(1,             # beta0 (intercept)
         0,             # beta1
         2)             # beta2
  },
  gen_ivs = function(n){
    X1 <- stats::rnorm(n, mean = 0, sd = 15)
    X2 <- stats::rnorm(n, 0, 3)

    #setup heteroskedastic errors
    grps = 5
    base_sd = 10
    sd_step = 5
    sequences <- split(1:n, ceiling((1:n)/(n/grps)))
    epsilon <- vector(mode = "numeric", length = n)

    # Generate increasing SDs for each group
    sds <- seq(from = base_sd, by = sd_step, length.out = grps)

    # Loop through groups to assign errors
    for(i in 1:grps) {
      group_X1 <- abs(X1[sequences[[i]]])
      epsilon[sequences[[i]]] <- stats::rnorm(n/grps, 0, group_X1*sds[i])
    }
    return(list(X1 = X1, X2 = X2, epsilon = epsilon))
  },
  b = 100,
  m = 10,
  mpat = function(mdat){ # define missing data pattern
    n <- nrow(mdat) * ncol(mdat)
    c1 <- stats::rbinom(n, 1, 0.7)
    c2 <- stats::rbinom(n, 1, 0.7)
    pattern_vec <- c1 * c2
    matrix(pattern_vec, nrow = nrow(mdat), ncol = ncol(mdat))
  },
  remove.collinear=FALSE,
  combine = mean
)

usethis::use_data(type1_2var_uncor_n75_het, overwrite = TRUE)
