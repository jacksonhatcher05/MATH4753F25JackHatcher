library(testthat)
library(MATH4753F25JackHatcher)

test_that("myncurve returns correct names", {
  result <- myncurve(0,1,1)
  expect_true(all(c("mu","sigma","probability") %in% names(result)))
})

test_that("probability is correct", {
  result <- myncurve(0,1,1)
  expect_equal(round(result$probability, 4), round(pnorm(1),4))
})
