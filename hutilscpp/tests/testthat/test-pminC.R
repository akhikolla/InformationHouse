test_that("pminC int", {
  x <- c(1L, -3L, 2L, .Machine$integer.max)
  expect_identical(pminC(x, 3L), pmin.int(x, 3L))
  y <- c(x, 4L)
  res <- pmin.int(y, 0L)
  expect_identical(pminC(y, 0L, in_place = TRUE), res)
  expect_equal(pminC(y, 0L), y)
  expect_equal(pminC(y, 0), y)


})

test_that("pminC in-place", {
  o <- c(-1, 0, 1)
  expect_equal(pminC(o, 0.5), c(-1, 0, 0.5))
  expect_equal(o, c(-1, 0, 1))
  pminC(o, 0.5, in_place = TRUE)
  expect_equal(o, c(-1, 0, 0.5))


  oi <- -1:5
  pmin2 <- pmin(-1:5, 2L)
  expect_equal(pminC(oi, 2L), pmin(oi, 2L))
  pminC(oi, 2L, in_place = TRUE)
  expect_equal(oi, pmin2)
})


test_that("pminC Error handling", {
  expect_error(pminC(1:5, 0L, in_place = "not logical"))
  expect_error(pminC("abc", "abc"), regexp = "numeric")
  expect_error(pminC(0:1, "abc"), regexp = "numeric")
})

test_that("pminC_real_real", {
  x <- runif(100)
  y <- x[2]
  expect_equal(pminC(x, y), pmin(x, y))
  x <- c(NaN, x)
  expect_equal(pminC(x, y, keep_nas = TRUE), pmin(x, y))
})


