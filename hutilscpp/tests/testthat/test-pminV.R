test_that("pminV works", {
  y <- sample.int(1e6)
  x <- sample.int(1e6)
  expect_identical(pminV(x, y), base::pmin(x, y))

  xd <- as.double(x)
  yd <- as.double(y)
  expect_identical(pminV(xd, yd), base::pmin(xd, yd))
})

test_that("Error handling", {
  expect_error(pminV(1:5, 1:6), "same length")
  expect_error(pminV(list(1:2), 1), "list")
  expect_error(pminV(1, list(1:2)), "list")
  expect_error(pminV("1", 1), "numeric")
  expect_error(pminV(1, "1"), "numeric")
})

test_that("in_place = TRUE", {
  x <- copy(1:10)
  pminV(x, x - 1L, in_place = TRUE)
  expect_equal(x, 0:9)
})

test_that("swap_xy", {
  expect_equal(pminV(1:5, 2:6 + 0, dbl_ok = TRUE),
               pminV(1:5 + 0, 2:6, dbl_ok = TRUE))
})

