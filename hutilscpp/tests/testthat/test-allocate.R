test_that("allocate works", {
  o <- allocate0_int(10)
  expect_identical(o, integer(10))
  p <- allocate0_dbl(10)
  expect_identical(p, double(10))

  xc <- allocate0_except(10, 2, -2)
  o <- integer(10)
  o[3] <- -2L
  expect_identical(xc, o)



})
