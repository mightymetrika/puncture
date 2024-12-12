test_that("puncture_lm_sim works", {
  sim_results <- puncture_lm_sim(
    n = 100,
    sim_iter = 2,
    beta_gen = function() {
      list(1,             # beta0 (intercept)
           1,             # beta1
           3)             # beta2
    },
    gen_ivs = function(n) {
      X1 <- stats::rnorm(n, 0, 15)
      X2 <- 0.6 * X1 + sqrt(1 - 0.6^2) * stats::rnorm(n)
      epsilon <- stats::rnorm(n, 0, abs(X1))
      list(X1 = X1, X2 = X2, epsilon = epsilon)
    },
    b = 5,
    m = 3,
    remove.collinear=FALSE,
    combine = mean
  ) |> suppressWarnings()

  expect_equal(length(sim_results), 2)
})


