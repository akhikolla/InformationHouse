library(testthat)
context("Testing rho returns precision and recall")

test_that("testing rho returns", {
  rho_set <- matrix(c(rep(1, times = 80), rep(0, times=100)), byrow=T, ncol=2)
  rho_set_1 <- rho(rho_set)
  expect_equal(round(rho_set_1$recall, 2), 1)
  expect_equal(round(rho_set_1$precision, 2), 1)
  
  rho_ct <- matrix(c(30, 0, 0, 60), ncol=2)
  rho_ct_1 <- rho(rho_ct)
  expect_equal(round(rho_ct_1$recall, 2), 1)
  expect_equal(round(rho_ct_1$precision, 2), 1)
})


test_that("testing precision/recall values", {
  rho_ct <- matrix(c(5,2,3,3), ncol=2)
  rho_ct_1 <- rho(rho_ct)
  
  expect_equal(rho_ct_1$recall, (5 / (5 + 3)))
  expect_equal(rho_ct_1$precision, (5 / (5 + 2)))
})