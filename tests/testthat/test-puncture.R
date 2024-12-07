test_that("puncture works", {
  df <- cars
  names(df) <- c("y", "x")
  res <- puncture(df, k = 4)
  expect_equal(nrow(res), 4)
  expect_equal(ncol(res), 3)

})
