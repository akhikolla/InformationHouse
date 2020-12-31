context("Proximity generation")
library(gretel)

test_that("Mode 'ogpv' works", {
  # check for node costs and no node costs
  # expect_equal_to_reference(generate_proximities(BuchDarrah19, mode = "ogpv", p = 3, node_costs = c(2,3,2,1,0)), file = "tests/testthat/ogpv_prox_cache")
})

test_that("Mode 'oppv' works", {
  # check for odds scale = 2 and odds_scale by node
  # expect_equal_to_reference(generate_proximities(BuchDarrah19, mode = "oppv", odds_scale_by_node = c(1,2,3,2,1)), file = "tests/testthat/oppv_prox_cache")
})

test_that("Mode 'sconduct' works", {
  # expect_equal_to_reference(generate_proximities(BuchDarrah19, mode = "sconductivity"), file = "tests/testthat/scond_prox_cache")
})