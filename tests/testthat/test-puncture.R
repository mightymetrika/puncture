test_that("puncture works", {
  res <- puncture(cars, b = 4, form = "speed ~ dist", term = "dist",
                  stacks = 2, thin = 1,
                  remove.collinear=FALSE)
  expect_equal(nrow(res), 4)
  expect_equal(ncol(res), 4)

})
