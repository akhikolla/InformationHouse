library(testthat)
context("Testing getHandSetIndices")

test_that("verify lack of positives", {
  set <- matrix(c(rep(1, times = 10), rep(0, times=50)), byrow=T, ncol=2)
  br <- 0.3
  hsl <- 20
  expect_error(getHandSetIndices(set, hsl, br), regexp = "Not enough positives")
})

test_that("ignore positives", {
  set <- matrix(c(rep(0, times=50)), byrow=T, ncol=2)
  br <- 0
  hsl <- 20
  
  hs <- getHandSetIndices(set, hsl, br)
  expect_equal(length(hs), hsl)
})

test_that("verify baserate", {
  set <- matrix(c(rep(1, times = 20), rep(0, times=100)), byrow=T, ncol=2)
  br <- 0.3
  hsl <- 20
  hs_inds <- getHandSetIndices(set, hsl, br)
  hs <- set[hs_inds,]
  
  expected_positive = ceiling(hsl * br)
  expect_gte(sum(hs[,1]), expected_positive)
})

test_that("verify reordering of positives", {
  set <- matrix(c(rep(1, times = 20), rep(0, times=100)), byrow=T, ncol=2)
  br <- 0.2
  hsl <- 40
  hs_inds <- getHandSetIndices(set, hsl, br)
  hs <- set[hs_inds,]
  hs_rater_one_pos <- which(hs[,1] == 1)
  
  expected_positive = ceiling(hsl * br)
  expect_gte(length(hs_rater_one_pos), expected_positive)
  expect_true(any(hs_rater_one_pos > expected_positive))
})