library(testthat)
context("Testing rho")

test_that("testing rho errors", {
  expect_error(rho(0.90, OcSBaserate = NULL), regexp = "Must give a baserate when a kappa")
  expect_error(rho(0.90, OcSBaserate = 0.2, testSetLength = NULL), regexp = "Must give a testSetLength when a kappa")
})

test_that("testing perfect rhos", {
  perfect_rho_set <- matrix(c(rep(1, times = 80), rep(0, times=100)), byrow=T, ncol=2)
  perfect_rho <- rho(perfect_rho_set)
  expect_equal(round(perfect_rho$rho, 2), 0.00)
  
  perfect_rho_ct <- matrix(c(30, 0, 0, 60), ncol=2)
  perfect_rho_2 <- rho(perfect_rho_ct)
  expect_equal(round(perfect_rho_2$rho, 2), 0.00)
})

test_that("testing rhoCT errors", {
  error_ct <- matrix(c(-1, 0, 0, 50), ncol=2)
  expect_error(rhoCT(error_ct), regexp = "Values in Contingency Table must be positive")
})

test_that("testing rhoSet w custom baserate", {
  perfect_rho_set <- matrix(c(rep(1, times = 40), rep(0, times=100)), byrow=T, ncol=2)
  
  rho_br <- rhoSet(perfect_rho_set, OcSBaserate = 0)
  expect_equal(round(rho_br$rho), 1)
})
