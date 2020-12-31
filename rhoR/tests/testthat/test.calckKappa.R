library(testthat)
context("Test kappa calculations")

test_that("kappa is correct", {
  expect_equal(floor(kappa(contingencyTable) * 10000), 6250)
  expect_equal(floor(kappa(codeSet) * 10000), 6250)
  expect_equal(floor(calcKappa(contingencyTable, isSet = FALSE) * 10000), 6250)
})

test_that("kappa short-circuit", {
  full_agmt <- matrix(c(20, 0, 0, 0), ncol = 2)
  expect_equal(calcKappa(full_agmt, isSet = FALSE), 1)
  
  full_agmt <- matrix(c(0, 0, 0, 20), ncol = 2)
  expect_equal(calcKappa(full_agmt, isSet = FALSE), 1)
  
  neg_agmt <- matrix(c(-1, 0, 0, 20), ncol = 2)
  expect_error(kappaCT(neg_agmt), regexp = "Values in Contingency Table must be positive")
  
  bad_ct <- matrix(c(10, 0, 0, 20), ncol = 1)
  expect_error(kappaCT(bad_ct), regexp = "Incorrect number of dimensions: Contingency table must be 2x2")
  
  bad_ct <- matrix(c(1.3, 0, 0, 20), ncol = 2)
  expect_error(kappaCT(bad_ct), regexp = "Contingency table values must be positive integers.")
})

test_that("kappa against threshold", {
  expect_true(calcKappa(contingencyTable, isSet = FALSE, kappaThreshold = 0.60)$above)
  expect_false(calcKappa(contingencyTable, isSet = FALSE, kappaThreshold = 0.70)$above)
})

test_that("CT baserate", {
  full_agmt <- matrix(c(20, 0, 0, 20), ncol = 2)
  br <- baserate(full_agmt)
  expect_equal(br$firstBaserate, br$secondBaserate)
  
  
  quarter_br <- matrix(c(20, 40, 0, 20), ncol = 2)
  br <- baserate(quarter_br)
  expect_false(br$firstBaserate == br$secondBaserate)
  expect_true(br$firstBaserate == 0.25)
  expect_true(br$secondBaserate == 0.75)
})

test_that("CT baserate error", {
  test_set <- matrix(c(20, 40, -1, 20), ncol = 2)
  expect_error(baserateCT(test_set))
})

test_that("Verify recall", {
  test_k <- rhoR:::calcKappa(contingencyTable, isSet = FALSE)
  test_r <- rhoR:::getR(test_k, 0.2, 0.8)
  test_r_c <- rhoR:::recall(test_k, 0.2, 0.8)
  expect_equal(round(test_r,   3), .606)
  expect_equal(round(test_r_c, 3), .606)
  expect_equal(round(test_r_c, 3), round(test_r, 3))
})
