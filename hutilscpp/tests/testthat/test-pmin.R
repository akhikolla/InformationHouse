context("test-pmin")

test_that("Error handling", {
  expect_error(pminV(1:5, 1:6),
               regexp = "length")
})




