library(testthat)
context("Test getHandset")

test_that("default getHandSet tests", {
  set_length <- 10
  base_rate <- 0.1
  
  test_set<- getHandSet(codeSet, set_length, base_rate, returnSet = TRUE)
  
  test_set_br <- baserate(test_set)
  
  expect_equal(colnames(test_set), c("first_rater", "second_rater"))
  expect_gte(test_set_br$firstBaserate, base_rate)
})

test_that("getHandSet tests without enought positives", {
  set_length <- 40
  base_rate <- 0.8
  
  expect_error(getHandSet(codeSet, set_length, base_rate, returnSet = TRUE))
})

test_that("default getHandSet return kappa", {
  set_length <- 40
  base_rate <- 0
  
  # pull the full example rhoR::codeSet, with expected K of 0.625
  hand_set_kappa <- getHandSet(codeSet, set_length, base_rate)
  
  expect_equal(hand_set_kappa, 0.625)
})

test_that("default getHandCT tests", {
  set_length <- 10
  base_rate <- 0.1
  
  test_set<- getHandCT(contingencyTable, set_length, base_rate, as_kappa = FALSE)
  test_set_br <- baserate(test_set)
  
  expect_gte(test_set_br$firstBaserate, base_rate)
  expect_equal(sum(test_set), set_length)
  
  testthat::expect_error(getHandCT(contingencyTable, set_length, 0.6, as_kappa = FALSE), regexp = "Not enough positives")
  
  ct <- matrix(c(300,0, 0, 6000), nrow = 2, byrow = T)
  ct_k <- getHandCT(ct, 200, 0.4, as_kappa = TRUE)
  expect_equal(ct_k, 1)
})

test_that("getHandSet returns proper class", {
  set_length <- 10
  base_rate <- 0.1
  
  test_set<- getHandSet(codeSet, set_length, base_rate, returnSet = TRUE)
  
  expect_is(test_set, "code.set")
})