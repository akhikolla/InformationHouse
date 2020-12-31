library(testthat)
library(triangulr)

################################################################################
## Setup

ptri_test <- function(q, a, b, c) {
  p <- (q - a)^2 / ((b - a) * (c - a))
  p[q <= a] <- 0
  p[b <= q] <- 1
  i <- c < q & q < b
  p[i] <- 1 - (b - q[i])^2 / ((b - a) * (b - c))
  p
}

ptri_vec_test <- function(q, a, b, c) {
  p <- (q - a)^2 / ((b - a) * (c - a))
  p[q <= a] <- 0
  p[b <= q] <- 1
  i <- c < q & q < b
  p[i] <- 1 - (b - q[i])^2 / ((b - a[i]) * (b - c))
  p
}

################################################################################
## Test cases for the cumulative distribution function

test_that("scalar q, scalar params, symmetric", {
  p <- ptri(0.5, min = 0, max = 1, mode = 0.5)
  p_test <- ptri_test(0.5, 0, 1, 0.5)
  expect_equal(p, p_test)
})

test_that("scalar q, scalar params, non-symmetric", {
  p <- ptri(0.5, min = 0, max = 1, mode = 0.8)
  p_test <- ptri_test(0.5, 0, 1, 0.8)
  expect_equal(p, p_test)
})

test_that("scalar q, scalar params, non-symmetric, upper_tail", {
  p <- ptri(0.4, min = 0, max = 1, mode = 0.8, lower_tail = FALSE)
  p_test <- ptri_test(0.4, 0, 1, 0.8)
  expect_equal(1 - p, p_test)
})

test_that("scalar q, scalar params, non-symmetric, log_p", {
  p <- ptri(0.5, min = 0, max = 1, mode = 0.8, log_p = TRUE)
  p_test <- ptri_test(0.5, 0, 1, 0.8)
  expect_equal(exp(p), p_test)
})

test_that("scalar q, scalar params, non-symmetric, upper_tail, log_p", {
  p <- ptri(0.4, min = 0, max = 1, mode = 0.8, lower_tail = FALSE, log_p = TRUE)
  p_test <- ptri_test(0.4, 0, 1, 0.8)
  expect_equal(1 - exp(p), p_test)
})

test_that("vector q, scalar params, symmetric", {
  p <- ptri(c(0.1, 0.6, 0.9), min = 0, max = 1, mode = 0.5)
  p_test <- ptri_test(c(0.1, 0.6, 0.9), 0, 1, 0.5)
  expect_equal(p, p_test)
})

test_that("vector q, scalar params, non-symmetric", {
  p <- ptri(c(0.1, 0.6, 0.9), min = 0, max = 1, mode = 0.8)
  p_test <- ptri_test(c(0.1, 0.6, 0.9), 0, 1, 0.8)
  expect_equal(p, p_test)
})

test_that("vector q, scalar params, non-symmetric, log_p", {
  p <- ptri(c(0.1, 0.6, 0.9), min = 0, max = 1, mode = 0.8, log_p = TRUE)
  p_test <- ptri_test(c(0.1, 0.6, 0.9), 0, 1, 0.8)
  expect_equal(exp(p), p_test)
})

test_that("vector q, vector params, symmetric", {
  p <- ptri(c(0.1, 0.6, 0.9), min = c(0, 0, 0), max = 2, mode = 1)
  p_test <- ptri_vec_test(c(0.1, 0.6, 0.9), c(0, 0, 0), 2, 1)
  expect_equal(p, p_test)
})

test_that("vector q, vector params, non-symmetric", {
  p <- ptri(c(0.1, 0.6, 0.9), min = c(-0.1, 0, 0.1), max = 3, mode = 2)
  p_test <- ptri_vec_test(c(0.1, 0.6, 0.9), c(-0.1, 0, 0.1), 3, 2)
  expect_equal(p, p_test)
})

test_that("vector q, vector params, upper_tail", {
  p <- ptri(c(0.1, 0.6, 0.9), min = c(0, 0.01, 0.02), max = 1, mode = 0.8,
            lower_tail = FALSE)
  p_test <- ptri_vec_test(c(0.1, 0.6, 0.9), c(0, 0.01, 0.02), 1, 0.8)
  expect_equal(1 - p, p_test)
})

test_that("vector q, vector params, log_p", {
  p <- ptri(c(0.1, 0.6, 0.9), min = c(0, 0.01, 0.02), max = 1, mode = 0.8,
            log_p = TRUE)
  p_test <- ptri_vec_test(c(0.1, 0.6, 0.9), c(0, 0.01, 0.02), 1, 0.8)
  expect_equal(exp(p), p_test)
})

test_that("vector q, vector params, upper_tail, log_p", {
  p <- ptri(c(0.1, 0.6, 0.9), min = c(0, 0.01, 0.02), max = 1, mode = 0.8, lower_tail = FALSE, log_p = TRUE)
  p_test <- ptri_vec_test(c(0.1, 0.6, 0.9), c(0, 0.01, 0.02), 1, 0.8)
  expect_equal(1 - exp(p), p_test)
})

test_that("Mode at bound, min == mode", {
  q <- c(0.1, 0.6, 0.9)
  p <- ptri(q, min = 0, max = 1, mode = 0)
  p_test <- ptri_test(c(0.1, 0.6, 0.9), 0, 1, 0)
  expect_equal(p, p_test)
  p <- ptri(q, min = c(0, 0, 0), max = 1, mode = 0)
  p_test <- ptri_vec_test(c(0.1, 0.6, 0.9), c(0, 0, 0), 1, 0)
  expect_equal(p, p_test)
})

test_that("Mode at bound, max == mode", {
  p <- ptri(c(0.1, 0.6, 0.9), min = 0, max = 1, mode = 1)
  p_test <- ptri_test(c(0.1, 0.6, 0.9), 0, 1, 1)
  expect_equal(p, p_test)
  p <- ptri(c(0.1, 0.6, 0.9), min = c(0, 0, 0), max = 1, mode = 1)
  p_test <- ptri_vec_test(c(0.1, 0.6, 0.9), c(0, 0, 0), 1, 1)
  expect_equal(p, p_test)
})

test_that("Zeros produced, q < min", {
  p <- ptri(c(0.1, 0.6, 0.9), min = 0.4, max = 1.4, mode = 0.9)
  p_test <- ptri_test(c(0.1, 0.6, 0.9), 0.4, 1.4, 0.9)
  expect_equal(p, p_test)
  p <- ptri(c(0.1, 0.6, 0.9), min = c(0.4, 0.4, 0.4), max = 1.4, mode = 0.9)
  p_test <- ptri_vec_test(c(0.1, 0.6, 0.9), c(0.4, 0.4, 0.4), 1.4, 0.9)
  expect_equal(p, p_test)
})

test_that("Ones produced, q > max", {
  p <- ptri(c(0.1, 0.6, 0.9), min = -0.5, max = 0.5, mode = 0)
  p_test <- ptri_test(c(0.1, 0.6, 0.9), -0.5, 0.5, 0)
  expect_equal(p, p_test)
  p <- ptri(c(0.1, 0.6, 0.9), min = c(-0.5, -0.5, -0.5), max = 0.5, mode = 0)
  p_test <- ptri_vec_test(c(0.1, 0.6, 0.9), c(-0.5, -0.5, -0.5), 0.5, 0)
  expect_equal(p, p_test)
})

test_that("NaN produced, mode < min", {
  p <- expect_warning(ptri(1, min = 0, max = 2, mode = -1))
  expect_equal(p, NaN)
  p <- expect_warning(ptri(c(1, 2, 3), min = c(-1, 0, 1), max = 3, mode = -1))
  expect_equal(p, c(0.75, NaN, NaN))
})

test_that("NaN produced, min == mode == max", {
  p <- expect_warning(ptri(3, min = 3, max = 3, mode = 3))
  expect_equal(p, NaN)
  p <- expect_warning(ptri(c(1, 2, 3), min = c(0, 1, 3), max = 3, mode = 3))
  expect_equal(round(p, 7), c(0.1111111, 0.25, NaN))
})

test_that("NaN produced, min > max", {
  p <- expect_warning(ptri(1, min = 0, max = -1, mode = 1))
  expect_equal(p, NaN)
  p <- expect_warning(ptri(c(1, 2, 3), min = c(1, 2, 3), max = 2.5, mode = 2))
  expect_equal(p, c(0, 0, NaN))
})

test_that("Error, NULL arguments", {
  expect_error(ptri(q = NULL))
  expect_error(ptri(q = 1, min = NULL))
  expect_error(ptri(q = 1, max = NULL))
  expect_error(ptri(q = 1, mode = NULL))
})

test_that("Error, Non-numeric arguments", {
  expect_error(ptri(q = "1"))
  expect_error(ptri(q = 1, min = "0"))
  expect_error(ptri(q = 1, max = "1"))
  expect_error(ptri(q = 1, mode = "0.5"))
})

test_that("Error, Non-logical argument", {
  expect_error(ptri(q = 1, lower_tail = "TRUE"))
  expect_error(ptri(q = 1, log_p = "FALSE"))
})

test_that("Error, illegal recycling", {
  expect_error(ptri(seq(0.1, 1, 0.1), min = 0, max = c(1, 2), mode = 0.5))
  expect_error(ptri(seq(0.1, 1, 0.1), min = c(0, 0.1), max = 1, mode = 0.5))
  expect_error(ptri(seq(0.1, 1, 0.1), min = 0, max = 1, mode = c(0.5, 0.6)))
})
