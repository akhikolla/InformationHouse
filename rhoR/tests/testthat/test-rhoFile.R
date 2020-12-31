library(testthat)
context("Test rho file")

test_that("Reading ratings from file", {
  fileName <- system.file("extdata", "codeSet-multiple.csv", package = "rhoR", mustWork = TRUE)
  r_1 <- rho.file(fileName, col1 = "rater4", col2 = "rater3", ScSKappaThreshold = 0.6)
  r_2 <- rho.file(fileName, col1 = "rater1", col2 = "rater3", ScSKappaThreshold = 0.6)

  testthat::expect_equivalent(r_1$rho, 0.00, tolerance = 0.01)
  testthat::expect_equivalent(r_1$kappa, 1.00, tolerance = 0.01)
  testthat::expect_equivalent(r_2$kappa, 0.85, tolerance = 0.01)
  
  testthat::expect_error(rho.file(fileName, col1 = "rater1", col2 = "rater2"))
})
