library(testthat)
library(triangulr)

################################################################################
## Setup

qtri_test <- function(p, a, b, c) {
  q <- a + sqrt(p * (b - a) * (c - a))
  i <- p >= (c - a) / (b - a)
  q[i] <- b - sqrt((1 - p[i]) * (b - a) * (b - c))
  q
}

qtri_vec_test <- function(p, a, b, c) {
  q <- a + sqrt(p * (b - a) * (c - a))
  i <- p >= (c - a) / (b - a)
  q[i] <- b - sqrt((1 - p[i]) * (b - a[i]) * (b - c))
  q
}

################################################################################
## Test cases for the quantile function

test_that("scalar p, scalar params, symmetric", {
  q <- qtri(0.5, min = 0, max = 1, mode = 0.5)
  q_test <- qtri_test(0.5, 0, 1, 0.5)
  expect_equal(q, q_test)
})

test_that("scalar p, scalar params, non-symmetric", {
  q <- qtri(0.5, min = 0, max = 1, mode = 0.8)
  q_test <- qtri_test(0.5, 0, 1, 0.8)
  expect_equal(q, q_test)
})

test_that("scalar p, scalar params, non-symmetric, upper_tail", {
  q <- qtri(0.4, min = 0, max = 1, mode = 0.8, lower_tail = FALSE)
  q_test <- qtri_test(1 - 0.4, 0, 1, 0.8)
  expect_equal(q, q_test)
})

test_that("scalar p, scalar params, non-symmetric, log_p", {
  q <- qtri(log(0.4), min = 0, max = 1, mode = 0.8, log_p = TRUE)
  q_test <- qtri_test(0.4, 0, 1, 0.8)
  expect_equal(q, q_test)
})

test_that("scalar p, scalar params, non-symmetric, upper_tail, log_p", {
  q <- qtri(log(0.4), min = 0, max = 1, mode = 0.8, lower_tail = FALSE,
            log_p = TRUE)
  q_test <- qtri_test(1 - 0.4, 0, 1, 0.8)
  expect_equal(q, q_test)
})

test_that("vector p, scalar params, symmetric", {
  q <- qtri(c(0.1, 0.6, 0.9), min = 0, max = 1, mode = 0.5)
  q_test <- qtri_test(c(0.1, 0.6, 0.9), 0, 1, 0.5)
  expect_equal(q, q_test)
})

test_that("vector p, scalar params, non-symmetric", {
  q <- qtri(c(0.1, 0.6, 0.9), min = 0, max = 1, mode = 0.8)
  q_test <- qtri_test(c(0.1, 0.6, 0.9), 0, 1, 0.8)
  expect_equal(q, q_test)
})

test_that("vector p, scalar params, non-symmetric, log_p", {
  q <- qtri(log(c(0.1, 0.6, 0.9)), min = 0, max = 1, mode = 0.8, log_p = TRUE)
  q_test <- qtri_test(c(0.1, 0.6, 0.9), 0, 1, 0.8)
  expect_equal(q, q_test)
})

test_that("vector p, vector params", {
  q <- qtri(c(0.1, 0.6, 0.9), min = c(0, 0.05, 0.1), max = 3, mode = 2)
  q_test <- qtri_vec_test(c(0.1, 0.6, 0.9), c(0, 0.05, 0.1), 3, 2)
  expect_equal(q, q_test)
})

test_that("vector q, vector params, upper_tail", {
  q <- qtri(c(0.1, 0.6, 0.9), min = c(0, 0.05, 0.1), max = 3, mode = 2,
            lower_tail = FALSE)
  q_test <- qtri_vec_test(1 - c(0.1, 0.6, 0.9), c(0, 0.05, 0.1), 3, 2)
  expect_equal(q, q_test)
})

test_that("vector q, vector params, log_p", {
  q <- qtri(c(0.1, 0.6, 0.9), min = c(0, 0.05, 0.1), max = 3, mode = 2,
            log_p = FALSE)
  q_test <- qtri_vec_test(c(0.1, 0.6, 0.9), c(0, 0.05, 0.1), 3, 2)
  expect_equal(q, q_test)
})

test_that("vector q, vector params, upper_tail, log_p", {
  q <- qtri(log(c(0.1, 0.6, 0.9)), min = c(0, 0.05, 0.1), max = 3, mode = 2,
            lower_tail = FALSE, log_p = TRUE)
  q_test <- qtri_vec_test(1 - c(0.1, 0.6, 0.9), c(0, 0.05, 0.1), 3, 2)
  expect_equal(q, q_test)
})

test_that("Mode at bound, min == mode", {
  q <- qtri(c(0.1, 0.6, 0.9), min = 0, max = 1, mode = 0)
  q_test <- qtri_test(c(0.1, 0.6, 0.9), 0, 1, 0)
  expect_equal(q, q_test)
  q <- qtri(c(0.1, 0.6, 0.9), min = c(0, 0, 0), max = 1, mode = 0)
  q_test <- qtri_vec_test(c(0.1, 0.6, 0.9), c(0, 0, 0), 1, 0)
  expect_equal(q, q_test)
})

test_that("Mode at bound, max == mode", {
  q <- qtri(c(0.1, 0.6, 0.9), min = 0, max = 1, mode = 1)
  q_test <- qtri_test(c(0.1, 0.6, 0.9), 0, 1, 1)
  expect_equal(q, q_test)
  q <- qtri(c(0.1, 0.6, 0.9), min = c(0, 0.05, 0.1), max = 1, mode = 1)
  q_test <- qtri_vec_test(c(0.1, 0.6, 0.9), c(0, 0.05, 0.1), 1, 1)
  expect_equal(q, q_test)
})

test_that("NaN produced, p < 0 || p > 1", {
  q <- expect_warning(qtri(c(-0.5, 1.5), min = 0.4, max = 1.4, mode = 0.9))
  expect_equal(q, c(NaN, NaN))
  q <- expect_warning(
    qtri(c(-0.5, 1.5), min = c(0.4, 0.1), max = 1.4, mode = 0.9)
  )
  expect_equal(q, c(NaN, NaN))
})

test_that("NaN produced, mode < min", {
  q <- expect_warning(qtri(c(0.4, 0.5), min = 0, max = 4, mode = -1))
  expect_equal(q, c(NaN, NaN))
  q <- expect_warning(qtri(c(0.4, 0.5), min = 0, max = 1, mode = c(-1, 0.5)))
  expect_equal(q, c(NaN, 0.5))
})

test_that("NaN produced, min == mode == max", {
  q <- expect_warning(qtri(c(0.4, 0.5), min = 0, max = 0, mode = 0))
  expect_equal(q, c(NaN, NaN))
  q <- expect_warning(qtri(c(0, 0.5), min = 0, max = c(0, 1), mode = c(0, 0.5)))
  expect_equal(q, c(NaN, 0.5))
})

test_that("NaN produced, min > max", {
  q <- expect_warning(qtri(c(0.4, 0.5), min = 0.5, max = 0, mode = 0.5))
  expect_equal(q, c(NaN, NaN))
  q <- expect_warning(qtri(c(0, 0.5), min = c(1.5, 0), max = 1, mode = 0.5))
  expect_equal(q, c(NaN, 0.5))
})

test_that("Error, NULL arguments", {
  expect_error(qtri(p = NULL))
  expect_error(qtri(p = 1, min = NULL))
  expect_error(qtri(p = 1, max = NULL))
  expect_error(qtri(p = 1, mode = NULL))
})

test_that("Error, Non-numeric arguments", {
  expect_error(qtri(p = "1"))
  expect_error(qtri(p = 1, min = "0"))
  expect_error(qtri(p = 1, max = "1"))
  expect_error(qtri(p = 1, mode = "0.5"))
})

test_that("Error, Non-logical argument", {
  expect_error(qtri(p = 1, lower_tail = "TRUE"))
  expect_error(qtri(p = 1, log_p = "FALSE"))
})

test_that("Error, illegal recycling", {
  expect_error(qtri(c(0.1, 0.6, 0.9), min = 0, max = c(1, 2), mode = 0.5))
  expect_error(qtri(c(0.1, 0.6, 0.9), min = c(0, 0.1), max = 1, mode = 0.5))
  expect_error(qtri(c(0.1, 0.6, 0.9), min = 0, max = 1, mode = c(0.5, 0.6)))
})
