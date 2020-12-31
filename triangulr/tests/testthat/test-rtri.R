library(testthat)
library(triangulr)
library(dqrng)

################################################################################
## Setup

before_each <- function() set.seed(1)

rtri_rand <- function(n) {
  before_each()
  runif(n)
}

rtri_test <- function(n, a, b, c) {
  p <- rtri_rand(n)
  q <- a + sqrt(p * (b - a) * (c - a))
  i <- p >= (c - a) / (b - a)
  q[i] <- b - sqrt((1 - p[i]) * (b - a) * (b - c))
  q
}

rtri_vec_test <- function(n, a, b, c) {
  p <- rtri_rand(n)
  q <- a + sqrt(p * (b - a) * (c - a))
  i <- p >= (c - a) / (b - a)
  q[i] <- b - sqrt((1 - p[i]) * (b - a[i]) * (b - c))
  q
}

################################################################################
## Test cases for the random variate generator function

test_that("n == 1, scalar params, symmetric", {
  before_each()
  r <- rtri(1, min = 0, max = 1, mode = 0.5)
  r_test <- rtri_test(1, 0, 1, 0.5)
  expect_equal(r, r_test)
})

test_that("n == 1, scalar params, non-symmetric", {
  before_each()
  r <- rtri(1, min = 0, max = 1, mode = 0.8)
  r_test <- rtri_test(1, 0, 1, 0.8)
  expect_equal(r, r_test)
})

test_that("n > 1, scalar params, symmetric", {
  before_each()
  r <- rtri(3, min = 0, max = 1, mode = 0.5)
  r_test <- rtri_test(3, 0, 1, 0.5)
  expect_equal(r, r_test)
})

test_that("n > 1, scalar params, non-symmetric", {
  before_each()
  r <- rtri(3, min = 0, max = 1, mode = 0.8)
  r_test <- rtri_test(3, 0, 1, 0.8)
  expect_equal(r, r_test)
})

test_that("n > 1, vector params, symmetric", {
  before_each()
  r <- rtri(3, min = c(0, 1, 2), max = 3, mode = 2)
  r_test <- rtri_vec_test(3, c(0, 1, 2), 3, 2)
  expect_equal(r, r_test)
})

test_that("n > 1, vector params, non-symmetric", {
  before_each()
  r <- rtri(3, min = c(0.1, 0.6, 0.9), max = 4, mode = 3.8)
  r_test <- rtri_vec_test(3, c(0.1, 0.6, 0.9), 4, 3.8)
  expect_equal(r, r_test)
})

test_that("Mode at bound, min == mode", {
  before_each()
  r <- rtri(3, min = 0, max = 1, mode = 0)
  r_test <- rtri_test(3, 0, 1, 0)
  expect_equal(r, r_test)
  before_each()
  r <- rtri(3, min = c(0, 0, 0), max = 1, mode = 0)
  r_test <- rtri_vec_test(3, c(0, 0, 0), 1, 0)
  expect_equal(r, r_test)
})

test_that("Mode at bound, max == mode", {
  before_each()
  r <- rtri(3, min = 0, max = 1, mode = 1)
  r_test <- rtri_test(3, 0, 1, 1)
  expect_equal(r, r_test)
  before_each()
  r <- rtri(3, min = c(0.1, 0.5, 0.9), max = 1, mode = 1)
  r_test <- rtri_vec_test(3, c(0.1, 0.5, 0.9), 1, 1)
  expect_equal(r, r_test)
})

test_that("NaN produced, mode < min", {
  before_each()
  r <- expect_warning(rtri(2, min = 1, max = 2, mode = 0))
  expect_equal(r, c(NaN, NaN))
  before_each()
  r <- expect_warning(rtri(2, min = c(1, 1), max = 2, mode = 0))
  expect_equal(round(r, 7), c(NaN, NaN))
})

test_that("NaN produced, min == mode == max", {
  before_each()
  r <- expect_warning(rtri(2, min = 0, max = 0, mode = 0))
  expect_equal(r, c(NaN, NaN))
  before_each()
  r <- expect_warning(rtri(2, min = 0, max = c(0, 0), mode = 0))
  expect_equal(round(r, 7), c(NaN, NaN))
})

test_that("NaN produced, min > max", {
  before_each()
  r <- expect_warning(rtri(2, min = 0, max = -1, mode = 1))
  expect_equal(r, c(NaN, NaN))
  before_each()
  r <- expect_warning(rtri(2, min = 0, max = c(-1, -1), mode = 1))
  expect_equal(round(r, 7), c(NaN, NaN))
})

test_that("Error, Negative n", {
  before_each()
  expect_error(rtri(n = -1))
})

test_that("Error, NULL arguments", {
  before_each()
  expect_error(rtri(n = NULL))
  expect_error(rtri(n = 1, min = NULL))
  expect_error(rtri(n = 1, max = NULL))
  expect_error(rtri(n = 1, mode = NULL))
})

test_that("Error, Non-numeric arguments", {
  before_each()
  expect_error(rtri(n = "1"))
  expect_error(rtri(n = 1, min = "0"))
  expect_error(rtri(n = 1, max = "1"))
  expect_error(rtri(n = 1, mode = "0.5"))
})

test_that("Error, illegal recycling", {
  before_each()
  expect_error(rtri(n = 10, min = 0, max = c(1, 2), mode = 0.5))
  expect_error(rtri(n = 10, min = c(0, 0.1), max = 1, mode = 0.5))
  expect_error(rtri(n = 10, min = 0, max = 1, mode = c(0.5, 0.6)))
})

test_that("dqrunif, scalar params", {
  dqset.seed(1)
  r <- rtri(3, min = 0, max = 1, mode = 0.5, dqrng = TRUE)
  dqset.seed(1)
  p <- dqrunif(3)
  a <- 0
  b <- 1
  c <- 0.5
  q <- a + sqrt(p * (b - a) * (c - a))
  i <- p >= (c - a) / (b - a)
  q[i] <- b - sqrt((1 - p[i]) * (b - a) * (b - c))
  r_test <- q
  expect_equal(r, r_test)
})

test_that("dqrunif, vector params", {
  dqset.seed(1)
  r <- rtri(3, min = rep.int(0, 3), max = 1, mode = 0.5, dqrng = TRUE)
  dqset.seed(1)
  p <- dqrunif(3)
  a <- rep.int(0, 3)
  b <- 1
  c <- 0.5
  q <- a + sqrt(p * (b - a) * (c - a))
  i <- p >= (c - a) / (b - a)
  q[i] <- b - sqrt((1 - p[i]) * (b - a[i]) * (b - c))
  r_test <- q
  expect_equal(r, r_test)
})
