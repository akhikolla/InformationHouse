context("test-pmax_pure_c")

test_that("Error handling", {
  expect_error(.Call("do_c_pmax", 1:5, 1, 2), regexp = "is not a real")
})

test_that("Pure C works", {
  x <- c(-1, 0, 0.5)
  expect_equal(pmax_pure_c(x, 0), pmax(x, 0))
  expect_equal(pmax_pure_c(x, 0L), pmax(x, 0L))
  res <- .Call("do_c_pmax", x, -1, 1)
  expect_equal(res, x)
  res <- .Call("do_c_pmax", x, -0.25, 0.25)
  expect_equal(min(res), -0.25)
  expect_equal(max(res), 0.25)

  res <- .Call("do_c_pmax", x, Inf, 0.1)
  expect_equal(max(res), 0.1)
})

test_that("pmax_pure_c int", {
  x <- -5:5
  expect_equal(pmax_pure_c(x, 0L), pmax(x, 0L))
})

