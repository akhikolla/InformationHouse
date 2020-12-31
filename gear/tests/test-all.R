library(gear)
if (requireNamespace("testthat", quietly = TRUE)) {
  library(testthat)
  test_check("gear")
} else {
  message("tests could not be run because testthat is not installed")
}
