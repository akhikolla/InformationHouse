library(testthat)
library(triangulr)

################################################################################
## Setup

mgtri_test <- function(t, a, b, c) {
  numer <- (b - c) * exp(a * t) -
    (b - a) * exp(c * t) +
    (c - a) * exp(b * t)
  denom <- (b - a) * (c - a) * (b - c) * t^2
  2 * numer / denom
}

################################################################################
## Test cases for the moment generating function

test_that("scalar t, scalar params, symmetric", {
  mg <- mgtri(0.5, min = 0, max = 1, mode = 0.5)
  mg_test <- mgtri_test(0.5, 0, 1, 0.5)
  expect_equal(mg, mg_test)
})

test_that("scalar t, scalar params, non-symmetric", {
  mg <- mgtri(0.5, min = 0, max = 1, mode = 0.8)
  mg_test <- mgtri_test(0.5, 0, 1, 0.8)
  expect_equal(mg, mg_test)
})

test_that("vector t, scalar params, symmetric", {
  mg <- mgtri(1:3, min = 0, max = 1, mode = 0.5)
  mg_test <- mgtri_test(1:3, 0, 1, 0.5)
  expect_equal(mg, mg_test)
})

test_that("vector t, scalar params, non-symmetric", {
  mg <- mgtri(1:3, min = 0, max = 1, mode = 0.8)
  mg_test <- mgtri_test(1:3, 0, 1, 0.8)
  expect_equal(mg, mg_test)
})

test_that("vector t, vector params, symmetric", {
  mg <- mgtri(1:3, min = 0:2, max = 2:4, mode = 1:3)
  mg_test <- mgtri_test(1:3, 0:2, 2:4, 1:3)
  expect_equal(mg, mg_test)
})

test_that("vector t, vector params, non-symmetric", {
  mg <- mgtri(1:3, min = 0:2, max = 8:10, mode = 7:9)
  mg_test <- mgtri_test(1:3, 0:2, 8:10, 7:9)
  expect_equal(mg, mg_test)
})

test_that("vector t, vector params recycled, symmetric", {
  mg <- mgtri(1:3, min = c(0, 0, 0), max = 2, mode = 1)
  mg_test <- mgtri_test(1:3, c(0, 0, 0), 2, 1)
  expect_equal(mg, mg_test)
})

test_that("vector t, vector params recycled, non-symmetric", {
  mg <- mgtri(1:3, min = c(0, 0.5, 1), max = 3, mode = 2)
  mg_test <- mgtri_test(1:3, c(0, 0.5, 1), 3, 2)
  expect_equal(mg, mg_test)
})

test_that("NaN produced, t == 0", {
  mg <- expect_warning(mgtri(0:2, min = 0.4, max = 1.4, mode = 0.9))
  mg_test <- mgtri_test(0:2, 0.4, 1.4, 0.9)
  expect_equal(mg, mg_test)
  mg <- expect_warning(mgtri(0:2, min = c(0.4, 0.5, 0.6), max = 1.4, mode = 0.9))
  mg_test <- mgtri_test(0:2, c(0.4, 0.5, 0.6), 1.4, 0.9)
  expect_equal(mg, mg_test)
})

test_that("NaN produced, Mode at bound, min == mode", {
  mg <- expect_warning(mgtri(1:3, min = 0, max = 1, mode = 0))
  mg_test <- mgtri_test(1:3, 0, 1, 0)
  expect_equal(mg, mg_test)
  mg <- expect_warning(mgtri(1:3, min = 0:2, max = 3:5, mode = 2))
  mg_test <- mgtri_test(1:3, 0:2, 3:5, 2)
  expect_equal(mg, mg_test)
})

test_that("NaN produced, Mode at bound, max == mode", {
  mg <- expect_warning(mgtri(1:3, min = 0, max = 1, mode = 1))
  mg_test <- mgtri_test(1:3, 0, 1, 1)
  expect_equal(mg, mg_test)
  mg <- expect_warning(mgtri(1:3, min = 0:2, max = 3:5, mode = 3))
  mg_test <- mgtri_test(1:3, 0:2, 3:5, 3)
  expect_equal(mg, mg_test)
})

test_that("NaN produced, mode < min", {
  mg <- expect_warning(mgtri(1:3, min = 0, max = 2, mode = -1))
  expect_equal(mg, rep.int(NaN, 3))
  mg <- expect_warning(mgtri(1:3, min = c(0, 0.5, 1), max = 2, mode = -1))
  expect_equal(mg, rep.int(NaN, 3))
})

test_that("NaN produced, min == mode == max", {
  mg <- expect_warning(mgtri(1:3, min = 2, max = 2, mode = 2))
  expect_equal(mg, rep.int(NaN, 3))
  mg <- expect_warning(mgtri(1:3, min = c(0, 0.5, 1), max = 2, mode = 2))
  expect_equal(mg, rep.int(NaN, 3))
})

test_that("NaN produced, min > max", {
  mg <- expect_warning(mgtri(1:3, min = 0, max = -1, mode = 1))
  expect_equal(mg, rep.int(NaN, 3))
  mg <- expect_warning(mgtri(1:3, min = c(0, 0.5, 1), max = -1, mode = 1))
  expect_equal(mg, rep.int(NaN, 3))
})

test_that("Error, NULL arguments", {
  expect_error(mgtri(t = NULL))
  expect_error(mgtri(t = 1, min = NULL))
  expect_error(mgtri(t = 1, max = NULL))
  expect_error(mgtri(t = 1, mode = NULL))
})

test_that("Error, Non-numeric arguments", {
  expect_error(mgtri(t = "1"))
  expect_error(mgtri(t = 1, min = "0"))
  expect_error(mgtri(t = 1, max = "1"))
  expect_error(mgtri(t = 1, mode = "0.5"))
})

test_that("Error, illegal recycling", {
  expect_error(mgtri(seq(0.1, 1, 0.1), min = 0, max = c(1, 2), mode = 0.5))
  expect_error(mgtri(seq(0.1, 1, 0.1), min = c(0, 0.1), max = 1, mode = 0.5))
  expect_error(mgtri(seq(0.1, 1, 0.1), min = 0, max = 1, mode = c(0.5, 0.6)))
})
