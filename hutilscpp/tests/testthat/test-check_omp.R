test_that("check_omp works", {
  expect_error(check_omp(NA), "logical")
  expect_error(check_omp(NA_integer_), "nThread")
  expect_error(check_omp(5.5), "whole number")
  skip_if_not_installed("parallel")
  skip_if_not(has_openmp())
  expect_error(check_omp(parallel::detectCores() + 1L), "nThread")
})
