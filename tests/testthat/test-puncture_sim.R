# test_that("puncture_sim works", {
#   sim_results <- puncture_sim(
#     n = 75,
#     sim_iter = 5,
#     beta_gen = function() {
#       list(1,             # beta0 (intercept)
#            0,             # beta1
#            3)             # beta2
#     },
#     gen_ivs = function(n) {
#       X1 <- stats::rnorm(n, 0, 20)
#       X2 <- 0.6 * X1 + sqrt(1 - 0.6^2) * stats::rnorm(n)
#       epsilon <- stats::rnorm(n, 0, abs(X1))
#       list(X1 = X1, X2 = X2, epsilon = epsilon)
#     },
#     k = 50,
#     b = 10
#   )
#   expect_equal(length(sim_results), 2)
# })
