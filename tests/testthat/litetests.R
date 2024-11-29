library(testthat)
library(XICOR)

test_that("xicor computes correct value for known case", {
  set.seed(123)
  x <- rnorm(100)
  y <- x^2 + rnorm(100, sd = 0.1)
  result <- xicor(x, y)
  expect_equal(result$estimate, 0.5, tolerance = 0.1)
})

test_that("xicor handles NA values correctly", {
  x <- c(1, 2, NA, 4, 5)
  y <- c(1, 4, 9, NA, 25)
  expect_warning(xicor(x, y))
})

test_that("xicor throws error for non-numeric input", {
  expect_error(xicor(letters, 1:26))
})
