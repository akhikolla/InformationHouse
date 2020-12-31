context("Optimal paths")
library(gretel)


test_that("'opt_gpv' finds the optimal path", {
  # without and with node costs
  expect_equal(opt_gpv(BuchDarrah19, p = 1, source = 1, target = 5), c(1,3,5))
  expect_equal(opt_gpv(BuchDarrah19, p = 1, node_costs = c(0,0,5,0,0), source = 1, target = 5), c(1,2,4,5))
})

test_that("'all_opt_gpv' finds optimals paths", {
  # without and with node costs
  # expect_equal_to_reference(all_opt_gpv(BuchDarrah19, p = 0, node_costs = c(0,0,5,0,0)), file = "tests/testthat/ogpv_cache")
})

test_that("'opt_ppv' finds the optimal path", {
  # without and with nodewise odds scales
  expect_equal(opt_ppv(BuchDarrah19, source = 1, target = 5), c(1,3,5))
  expect_equal(opt_ppv(BuchDarrah19, odds_scale_by_node = c(1,1,9,1,1), source = 1, target = 5), c(1,5))
})

test_that("'all_opt_ppv' finds optimals paths", {
  # without and with nodewise odds scales
  # expect_equal_to_reference(all_opt_ppv(BuchDarrah19, odds_scale_by_node = c(1,1,9,1,1)), file = "tests/testthat/oppv_cache")
})

test_that("unpack returns a valid path", {
  expect_equal(unpack(c(2,4,5,NA,4), source = 4, target = 1), c(4,2,1))
})

test_that("all four optimal path functions fail gracefully", {
  # two checks each, one if bad input for path, one if bad input for specialty param
  expect_error(opt_gpv(BuchDarrah19, source = 1.5, target = 4))
  expect_error(all_opt_gpv(BuchDarrah19, p = -1))
  expect_error(opt_ppv(BuchDarrah19, source = "a", target = "b"))
  expect_error(all_opt_ppv(BuchDarrah19, odds_scale = 0))
})

test_that("unpack fails gracefully", {
  # target should be 4
  expect_equal(unpack(c(NA,4,5,NA,4), source = 2, target = 1), NA)
})