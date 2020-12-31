context("test-units2km")

test_that("units2km works", {
  expect_equal(units2km("30km"), 30)
  expect_equal(units2km("30.0 km"), 30)
  expect_equal(units2km("305 m"), 0.305)
})
