library(testthat)
context("Testing createSimulatedCodeSet")

test_that("verify error checks", {
  len <- 100
  br <- 0.2
  kmin <- 0.6
  kmax <- 0.9
  pmin <- 0.1
  pmax <- 0.9
  tries <- 50
  
  expect_error(createSimulatedCodeSet(-len, br, kmin, kmax, pmin, pmax, tries), regexp = "length.*?positive")
  expect_error(createSimulatedCodeSet(len, -br, kmin, kmax, pmin, pmax, tries), regexp = "baserate.*?positive")
  expect_error(createSimulatedCodeSet(len, br, kmin - 1, kmax, pmin, pmax, tries), regexp = "kappaMin.*?0 and 1")
  expect_error(createSimulatedCodeSet(len, br, kmin + 1, kmax, pmin, pmax, tries), regexp = "kappaMin.*?0 and 1")
  expect_error(createSimulatedCodeSet(len, br, kmin, kmax - 1, pmin, pmax, tries), regexp = "kappaMax.*?0 and 1")
  expect_error(createSimulatedCodeSet(len, br, kmin, kmax + 1, pmin, pmax, tries), regexp = "kappaMax.*?0 and 1")
  expect_error(createSimulatedCodeSet(len, br, kappaMin = 0.9, kappaMax = 0.7, pmin, pmax, tries), regexp = "kappaMin.*?less.*?Max")
  expect_error(createSimulatedCodeSet(len, br, kmin, kmax, pmin - 1, pmax, tries), regexp = "precisionMin.*?0 and 1")
  expect_error(createSimulatedCodeSet(len, br, kmin, kmax, pmin + 1, pmax, tries), regexp = "precisionMin.*?0 and 1")
  expect_error(createSimulatedCodeSet(len, br, kmin, kmax, pmin, pmax - 1, tries), regexp = "precisionMax.*?0 and 1")
  expect_error(createSimulatedCodeSet(len, br, kmin, kmax, pmin, pmax + 1, tries), regexp = "precisionMax.*?0 and 1")
  expect_error(createSimulatedCodeSet(len, br, kmin, kmax, precisionMin = 0.9, precisionMax = 0.7, tries), regexp = "precisionMin.*less.*?Max")
  expect_error(createSimulatedCodeSet(len, br, kmin, kmax, pmin, pmax, tries = 0), regexp = "tries.*?1")
})

test_that("verify good set", {
  len <- 100
  br <- 0.2
  kmin <- 0.6
  kmax <- 0.9
  pmin <- 0.1
  pmax <- 0.9
  tries <- 50
  
  set <- createSimulatedCodeSet(len, br, kmin, kmax, pmin, pmax, tries)
  set_br <- baserate(set)
  set_k  <- kappa(set)
    
  expect_gte(set_br$firstBaserate, br)
  expect_gte(set_k, kmin)
  expect_lte(set_k, kmax)
})


test_that("impossible set", {
  len <- 100
  br <- 0.05
  kmin <- 0.10
  kmax <- 0.11
  pmin <- 0.6
  pmax <- 0.8
  tries <- 500
  
  expect_error(createSimulatedCodeSet(len, br, kmin, kmax, pmin, pmax, tries), regexp = "Unable to create.*26471")
})

test_that("explicit any_equal checks", {
  expect_false(any_equal(c(1.1001, 2.2), 11/10))
  expect_true(any_equal(c(1.1000000001, 2.2), 11/10))
  expect_true(any_equal(c(1.1, 2.2), 11/10))
  expect_true(any_equal(c(1.1, 2.2), 22/10))
})