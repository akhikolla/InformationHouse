library(testthat)
library(triangulr)

################################################################################
## Setup

ctri_test <- function(t, a, b, c) {
  n <- (b - c) * exp(1i * a * t) -
    (b - a) * exp(1i * c * t) +
    (c - a) * exp(1i * b * t)
  d <- (b - a) * (c - a) * (b - c) * t^2
  -2 * n / d
}

################################################################################
## Test cases for the characteristic function

test_that("scalar t, scalar params, symmetric", {
  c <- ctri(0.5, min = 0, max = 1, mode = 0.5)
  c_test <- ctri_test(0.5, 0, 1, 0.5)
  expect_equal(c, c_test)
})

test_that("scalar t, scalar params, non-symmetric", {
  c <- ctri(0.5, min = 0, max = 1, mode = 0.8)
  c_test <- ctri_test(0.5, 0, 1, 0.8)
  expect_equal(c, c_test)
})

test_that("vector t, scalar params, symmetric", {
  c <- ctri(1:3, min = 0, max = 1, mode = 0.5)
  c_test <- ctri_test(1:3, 0, 1, 0.5)
  expect_equal(c, c_test)
})

test_that("vector t, scalar params, non-symmetric", {
  c <- ctri(1:3, min = 0, max = 1, mode = 0.8)
  c_test <- ctri_test(1:3, 0, 1, 0.8)
  expect_equal(c, c_test)
})

test_that("vector t, vector params, symmetric", {
  c <- ctri(1:3, min = 0:2, max = 2:4, mode = 1:3)
  c_test <- ctri_test(1:3, 0:2, 2:4, 1:3)
  expect_equal(c, c_test)
})

test_that("vector t, vector params, non-symmetric", {
  c <- ctri(1:3, min = 0:2, max = 8:10, mode = 7:9)
  c_test <- ctri_test(1:3, 0:2, 8:10, 7:9)
  expect_equal(c, c_test)
})

test_that("vector t, vector params recycled, symmetric", {
  c <- ctri(1:3, min = c(0, 0, 0), max = 2, mode = 1)
  c_test <- ctri_test(1:3, c(0, 0, 0), 2, 1)
  expect_equal(c, c_test)
})

test_that("vector t, vector params recycled, non-symmetric", {
  c <- ctri(1:3, min = c(0, 0.5, 1), max = 3, mode = 2)
  c_test <- ctri_test(1:3, c(0, 0.5, 1), 3, 2)
  expect_equal(c, c_test)
})

test_that("NaN produced, t == 0", {
  c <- expect_warning(ctri(0:2, min = 0.4, max = 1.4, mode = 0.9))
  c_test <- ctri_test(0:2, 0.4, 1.4, 0.9)
  expect_equal(c, c_test)
  c <- expect_warning(ctri(0:2, min = c(0.4, 0.5, 0.6), max = 1.4, mode = 0.9))
  c_test <- ctri_test(0:2, c(0.4, 0.5, 0.6), 1.4, 0.9)
  expect_equal(c, c_test)
})

test_that("NaN produced, Mode at bound, min == mode", {
  c <- expect_warning(ctri(1:3, min = 0, max = 1, mode = 0))
  c_test <- ctri_test(1:3, 0, 1, 0)
  expect_equal(c, c_test)
  c <- expect_warning(ctri(1:3, min = 0:2, max = 3:5, mode = 2))
  c_test <- ctri_test(1:3, 0:2, 3:5, 2)
  expect_equal(c, c_test)
})

test_that("NaN produced, Mode at bound, max == mode", {
  c <- expect_warning(ctri(1:3, min = 0, max = 1, mode = 1))
  c_test <- ctri_test(1:3, 0, 1, 1)
  expect_equal(c, c_test)
  c <- expect_warning(ctri(1:3, min = 0:2, max = 3:5, mode = 3))
  c_test <- ctri_test(1:3, 0:2, 3:5, 3)
  expect_equal(c, c_test)
})

test_that("NaN produced, mode < min", {
  c <- expect_warning(ctri(1:3, min = 0, max = 2, mode = -1))
  expect_equal(c, rep.int(NaN + 0i, 3))
  c <- expect_warning(ctri(1:3, min = c(0, 0.5, 1), max = 2, mode = -1))
  expect_equal(c, rep.int(NaN + 0i, 3))
})

test_that("NaN produced, min == mode == max", {
  c <- expect_warning(ctri(1:3, min = 2, max = 2, mode = 2))
  expect_equal(c, rep.int(NaN + 0i, 3))
  c <- expect_warning(ctri(1:3, min = c(0, 0.5, 1), max = 2, mode = 2))
  expect_equal(c, rep.int(NaN + 0i, 3))
})

test_that("NaN produced, min > max", {
  c <- expect_warning(ctri(1:3, min = 0, max = -1, mode = 1))
  expect_equal(c, rep.int(NaN + 0i, 3))
  c <- expect_warning(ctri(1:3, min = c(0, 0.5, 1), max = -1, mode = 1))
  expect_equal(c, rep.int(NaN + 0i, 3))
})

test_that("Error, NULL arguments", {
  expect_error(ctri(t = NULL))
  expect_error(ctri(t = 1, min = NULL))
  expect_error(ctri(t = 1, max = NULL))
  expect_error(ctri(t = 1, mode = NULL))
})

test_that("Error, Non-numeric arguments", {
  expect_error(ctri(t = "1"))
  expect_error(ctri(t = 1, min = "0"))
  expect_error(ctri(t = 1, max = "1"))
  expect_error(ctri(t = 1, mode = "0.5"))
})

test_that("Error, illegal recycling", {
  expect_error(ctri(seq(0.1, 1, 0.1), min = 0, max = c(1, 2), mode = 0.5))
  expect_error(ctri(seq(0.1, 1, 0.1), min = c(0, 0.1), max = 1, mode = 0.5))
  expect_error(ctri(seq(0.1, 1, 0.1), min = 0, max = 1, mode = c(0.5, 0.6)))
})
