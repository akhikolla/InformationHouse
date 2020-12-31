library(testthat)
library(triangulr)

################################################################################
## Setup

# This function approximates the integral value
estri_test <- function(p, min, max, mode) {
  es <- function(x, a, b, c) {
    integrate(qtri, lower = 0, upper = x, min = a, max = b, mode = c)$value / x
  }
  mapply(es, x = p, a = min, b = max, c = mode)
}

expect_equal2 <- function(x, y, r = 5) expect_equal(round(x, r), round(y, r))

################################################################################
## Test cases for the expected shortfall function

test_that("scalar p, scalar params, symmetric", {
  es <- estri(0.5, min = 0, max = 1, mode = 0.5)
  es_test <- estri_test(0.5, 0, 1, 0.5)
  expect_equal2(es, es_test)
})

test_that("scalar p, scalar params, non-symmetric", {
  es <- estri(0.5, min = 0, max = 1, mode = 0.8)
  es_test <- estri_test(0.5, 0, 1, 0.8)
  expect_equal2(es, es_test)
})

test_that("scalar p, scalar params, non-symmetric, upper_tail", {
  es <- estri(0.4, min = 0, max = 1, mode = 0.8, lower_tail = FALSE)
  es_test <- estri_test(1 - 0.4, 0, 1, 0.8)
  expect_equal2(es, es_test)
})

test_that("scalar p, scalar params, non-symmetric, log_p", {
  es <- estri(log(0.5), min = 0, max = 1, mode = 0.8, log_p = TRUE)
  es_test <- estri_test(0.5, 0, 1, 0.8)
  expect_equal2(es, estri_test(0.5, 0, 1, 0.8))
})

test_that("scalar p, scalar params, non-symmetric, upper_tail, log_p", {
  es <- estri(log(0.4), min = 0, max = 1, mode = 0.8, lower_tail = FALSE,
              log_p = TRUE)
  es_test <- estri_test(1 - 0.4, 0, 1, 0.8)
  expect_equal2(es, es_test)
})

test_that("vector p, scalar params, symmetric", {
  es <- estri(c(0.1, 0.6, 0.9), min = 0, max = 1, mode = 0.5)
  es_test <- estri_test(c(0.1, 0.6, 0.9), 0, 1, 0.5)
  expect_equal2(es, es_test)
})

test_that("vector p, scalar params, non-symmetric", {
  es <- estri(c(0.1, 0.6, 0.9), min = 0, max = 1, mode = 0.8)
  es_test <- estri_test(c(0.1, 0.6, 0.9), 0, 1, 0.8)
  expect_equal2(es, es_test)
})

test_that("vector p, scalar params, non-symmetric, log_p", {
  p <- seq(0.1, 1, 0.1)
  es <- estri(log(c(0.1, 0.6, 0.9)), min = 0, max = 1, mode = 0.8, log_p = TRUE)
  es_test <- estri_test(c(0.1, 0.6, 0.9), 0, 1, 0.8)
  expect_equal2(es, es_test)
})

test_that("vector p, vector params, symmetric", {
  es <- estri(c(0.1, 0.6, 0.9), min = 0:2, max = 2:4, mode = 1:3)
  es_test <- estri_test(c(0.1, 0.6, 0.9), 0:2, 2:4, 1:3)
  expect_equal2(es, es_test)
})

test_that("vector p, vector params, non-symmetric", {
  es <- estri(c(0.1, 0.6, 0.9), min = 0:2, max = 5:7, mode = 4:6)
  es_test <- estri_test(c(0.1, 0.6, 0.9), 0:2, 5:7, 4:6)
  expect_equal2(es, es_test)
})

test_that("vector p, vector params, non-symmetric, upper_tail, log_p", {
  es <- estri(log(c(0.1, 0.6, 0.9)), min = 0:2, max = 5:7, mode = 4:6,
              lower_tail = FALSE, log_p = TRUE)
  es_test <- estri_test(1 - c(0.1, 0.6, 0.9), 0:2, 5:7, 4:6)
  expect_equal2(es, es_test)
})

test_that("vector p, vector params recycled, symmetric", {
  es <- estri(c(0.1, 0.6, 0.9), min = c(0, 0, 0), max = 1, mode = 0.5)
  es_test <- estri_test(c(0.1, 0.6, 0.9), c(0, 0, 0), 1, 0.5)
  expect_equal2(es, es_test)
})

test_that("vector p, vector params recycled, non-symmetric", {
  es <- estri(c(0.1, 0.6, 0.9), min = 1:3, max = 6:8, mode = 5:7)
  es_test <- estri_test(c(0.1, 0.6, 0.9), 1:3, 6:8, 5:7)
  expect_equal2(es, es_test)
})

test_that("Mode at bound, min == mode", {
  es <- estri(c(0.1, 0.6, 0.9), min = 0, max = 1, mode = 0)
  es_test <- estri_test(c(0.1, 0.6, 0.9), 0, 1, 0)
  expect_equal2(es, es_test)
  es <- estri(c(0.1, 0.6, 0.9), min = c(0, 0, 0), max = 1, mode = 0)
  es_test <- estri_test(c(0.1, 0.6, 0.9), 0, 1, 0)
  expect_equal2(es, es_test)
})

test_that("Mode at bound, max == mode", {
  es <- estri(c(0.1, 0.6, 0.9), min = 0, max = 1, mode = 1)
  es_test <- estri_test(c(0.1, 0.6, 0.9), 0, 1, 1)
  expect_equal2(es, es_test)
  es <- estri(c(0.1, 0.6, 0.9), min = c(0, 0, 0), max = 1, mode = 1)
  es_test <- estri_test(c(0.1, 0.6, 0.9), c(0, 0, 0), 1, 1)
  expect_equal2(es, es_test)
})

test_that("NaN produced, p == 0 || p < 0 || p > 1", {
  p <- c(-0.01, 0, 1.01)
  es <- expect_warning(estri(p, min = 0.4, max = 1.4, mode = 0.9))
  expect_equal(es, c(NaN, NaN, NaN))
  a <- c(0.3, 0.4, 0.5)
  es <- expect_warning(estri(p, min = a, max = 1.4, mode = 0.9))
  expect_equal(es, c(NaN, NaN, NaN))
})

test_that("NaN produced, mode < min", {
  es <- expect_warning(estri(0.5, min = 0, max = 2, mode = -1))
  expect_equal(es, NaN)
  es <- expect_warning(estri(c(0.5, 0.6), min = 0, max = 2, mode = c(-1, 1)))
  expect_equal2(es, c(NaN, 0.7308565))
})

test_that("NaN produced, min == mode == max", {
  es <- expect_warning(estri(0.5, min = 0.5, max = 0.5, mode = 0.5))
  expect_equal2(es, NaN)
  es <- expect_warning(
    estri(c(0.1, 0.5, 1), min = c(0, 1, 3), max = c(2, 3, 3), mode = c(1, 2, 3))
  )
  expect_equal2(es, c(0.2981424, 1.6666667, NaN))
})

test_that("NaN produced, min > max", {
  es <- expect_warning(estri(0.5, min = 1, max = 0, mode = 0))
  expect_equal(es, NaN)
  es <- expect_warning(
    estri(c(0, 0.5, 1), min = c(0, 1, 2), max = c(-1, 3, 4), mode = c(1, 2, 3))
  )
  expect_equal2(es, c(NaN, 1.666667, 3))
})

test_that("Error, NULL arguments", {
  expect_error(estri(p = NULL))
  expect_error(estri(p = 1, min = NULL))
  expect_error(estri(p = 1, max = NULL))
  expect_error(estri(p = 1, mode = NULL))
})

test_that("Error, Non-numeric arguments", {
  expect_error(estri(p = "1"))
  expect_error(estri(p = 1, min = "0"))
  expect_error(estri(p = 1, max = "1"))
  expect_error(estri(p = 1, mode = "0.5"))
})

test_that("Error, Non-logical argument", {
  expect_error(estri(p = 1, lower_tail = "TRUE"))
  expect_error(estri(p = 1, log_p = "FALSE"))
})

test_that("Error, illegal recycling", {
  p <- seq(0.1, 1, 0.1)
  expect_error(estri(p, min = 0, max = c(1, 2), mode = 0.5))
  expect_error(estri(p, min = c(0, 0.1), max = 1, mode = 0.5))
  expect_error(estri(p, min = 0, max = 1, mode = c(0.5, 0.6)))
})
