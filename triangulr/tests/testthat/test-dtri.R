library(testthat)
library(triangulr)

################################################################################
## Setup

dtri_test <- function(x, a, b, c) {
  p <- 2 * (x - a) / ((b - a) * (c - a))
  p[x == c] <- 2 / (b - a)
  m <- c < x & x <= b
  p[m] <- 2 * (b - x[m]) / ((b - a) * (b - c))
  p[x < a | b < x] <- 0
  p
}

dtri_vec_test <- function(x, a, b, c) {
  p <- 2 * (x - a) / ((b - a) * (c - a))
  i <- x == c
  p[i] <- 2 / (b - a[i])
  i <- c < x & x <= b
  p[i] <- 2 * (b - x[i]) / ((b - a[i]) * (b - c))
  p[x < a | b < x] <- 0
  p
}

################################################################################
## Test cases for the probability density function

test_that("scalar x, scalar params, symmetric", {
  d <- dtri(0.5, min = 0, max = 1, mode = 0.5)
  d_test <- dtri_test(0.5, 0, 1, 0.5)
  expect_equal(d, d_test)
})

test_that("scalar x, scalar params, non-symmetric", {
  d <- dtri(0.5, min = 0, max = 1, mode = 0.8)
  d_test <- dtri_test(0.5, 0, 1, 0.8)
  expect_equal(d, d_test)
})

test_that("scalar x, scalar params, non-symmetric, log-scale", {
  d <- dtri(0.5, min = 0, max = 1, mode = 0.8, log = TRUE)
  d_test <- dtri_test(0.5, 0, 1, 0.8)
  expect_equal(exp(d), d_test)
})

test_that("vector x, scalar params, symmetric", {
  d <- dtri(c(0.1, 0.6, 0.9), min = 0, max = 1, mode = 0.5)
  d_test <- dtri_test(c(0.1, 0.6, 0.9), 0, 1, 0.5)
  expect_equal(d, d_test)
})

test_that("vector x, scalar params, non-symmetric", {
  d <- dtri(c(0.1, 0.6, 0.9), min = 0, max = 1, mode = 0.8)
  d_test <- dtri_test(c(0.1, 0.6, 0.9), 0, 1, 0.8)
  expect_equal(d, d_test)
})

test_that("vector x, scalar params, non-symmetric, log-scale", {
  d <- dtri(c(0.1, 0.6, 0.9), min = 0, max = 1, mode = 0.8, log = TRUE)
  d_test <- dtri_test(c(0.1, 0.6, 0.9), 0, 1, 0.8)
  expect_equal(exp(d), d_test)
})

test_that("vector x, vector params, symmetric", {
  d <- dtri(c(0.1, 0.6, 0.9), min = c(0, 0, 0), max = 2, mode = 1)
  d_test <- dtri_vec_test(c(0.1, 0.6, 0.9), c(0, 0, 0), 2, 1)
  expect_equal(d, d_test)
})

test_that("vector x, vector params, non-symmetric", {
  d <- dtri(c(0.1, 0.6, 0.9), min = 0:2, max = 7, mode = 6)
  d_test <- dtri_vec_test(c(0.1, 0.6, 0.9), 0:2, 7, 6)
  expect_equal(d, d_test)
})

test_that("vector x, vector params, non-symmetric, log-scale", {
  d <- dtri(c(0.1, 0.6, 0.9), min = 0:2, max = 7, mode = 6, log = TRUE)
  d_test <- dtri_vec_test(c(0.1, 0.6, 0.9), 0:2, 7, 6)
  expect_equal(exp(d), d_test)
})

test_that("Mode at bound, min == mode", {
  d <- dtri(c(0.1, 0.6, 0.9), min = 0, max = 1, mode = 0)
  d_test <- dtri_test(c(0.1, 0.6, 0.9), 0, 1, 0)
  expect_equal(d, d_test)
  d <- dtri(c(0.1, 0.6, 0.9), min = 0:2, max = 3, mode = 2)
  d_test <- dtri_vec_test(c(0.1, 0.6, 0.9), 0:2, 3, 2)
  expect_equal(d, d_test)
})

test_that("Mode at bound, max == mode", {
  d <- dtri(c(0.1, 0.6, 0.9), min = 0, max = 1, mode = 1)
  d_test <- dtri_test(c(0.1, 0.6, 0.9), 0, 1, 1)
  expect_equal(d, d_test)
  d <- dtri(c(0.1, 0.6, 0.9), min = 0:2, max = 3, mode = 3)
  d_test <- dtri_vec_test(c(0.1, 0.6, 0.9), 0:2, 3, 3)
  expect_equal(d, d_test)
})

test_that("Zeros produced, x < min", {
  d <- dtri(1:3, min = 2, max = 4, mode = 3)
  d_test <- dtri_test(1:3, 2, 4, 3)
  expect_equal(d, d_test)
  d <- dtri(1:3, min = 2:4, max = 6, mode = 5)
  d_test <- dtri_vec_test(1:3, 2:4, 6, 5)
  expect_equal(d, d_test)
})

test_that("Zeros produced, x > max", {
  d <- dtri(1:3, min = 0, max = 2, mode = 1)
  d_test <- dtri_test(1:3, 0, 2, 1)
  expect_equal(d, d_test)
  d <- dtri(2:4, min = 0:2, max = 3, mode = 2)
  d_test <- dtri_vec_test(2:4, 0:2, 3, 2)
  expect_equal(d, d_test)
})

test_that("NaN produced, mode < min", {
  d <- expect_warning(dtri(1, min = 0, max = 2, mode = -1))
  expect_equal(d, NaN)
  d <- expect_warning(dtri(1:3, min = 1:3, max = 4, mode = 2))
  expect_equal(d, c(0, 1, NaN))
})

test_that("NaN produced, min == mode == max", {
  d <- expect_warning(dtri(1, min = 3, max = 3, mode = 3))
  expect_equal(d, NaN)
  d <- expect_warning(dtri(0:2, min = 0:2, max = 2, mode = 2))
  expect_equal(d, c(0, 0, NaN))
})

test_that("NaN produced, min > max", {
  d <- expect_warning(dtri(1, min = 0, max = -1, mode = 1))
  expect_equal(d, NaN)
  d <- expect_warning(dtri(c(1, 1.5, 2), min = c(1, 1, 3), max = 2, mode = 2))
  expect_equal(d, c(0, 1, NaN))
})

test_that("Error, NULL arguments", {
  expect_error(dtri(x = NULL))
  expect_error(dtri(x = 1, min = NULL))
  expect_error(dtri(x = 1, max = NULL))
  expect_error(dtri(x = 1, mode = NULL))
})

test_that("Error, Non-numeric arguments", {
  expect_error(dtri(x = "1"))
  expect_error(dtri(x = 1, min = "0"))
  expect_error(dtri(x = 1, max = "1"))
  expect_error(dtri(x = 1, mode = "0.5"))
})

test_that("Error, Non-logical argument", {
  expect_error(dtri(x = 1, log = "FALSE"))
})

test_that("Error, illegal recycling", {
  expect_error(dtri(seq(0.1, 1, 0.1), min = 0, max = c(1, 2), mode = 0.5))
  expect_error(dtri(seq(0.1, 1, 0.1), min = c(0, 0.1), max = 1, mode = 0.5))
  expect_error(dtri(seq(0.1, 1, 0.1), min = 0, max = 1, mode = c(0.5, 0.6)))
})
