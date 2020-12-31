context("test-range-cpp")

test_that("range works", {
  y <- 1
  expect_identical(range(y), range_rcpp(y)[1:2])
  y <- c(1, 1)
  expect_identical(range(y), range_rcpp(y)[1:2])
  y <- as.integer(y)
  expect_identical(range(y), range_rcpp(y)[1:2])
  y <- runif(500)
  expect_identical(range(y), range_rcpp(y)[1:2])
  y <- sample(1:500, size = 450)
  expect_identical(range(y), range_rcpp(y)[1:2])
  y <- seq(0, 100, length.out = 5)
  expect_identical(range(y), range_rcpp(y)[1:2])
  y <- rev(y)
  expect_identical(range(y), range_rcpp(y)[1:2])
  y <- as.integer(y)
  expect_identical(range(y), range_rcpp(y)[1:2])
  y <- rev(y)
  expect_identical(range(y), range_rcpp(y)[1:2])
})

test_that("range falls through", {
  y <- as.character(c(letters, 1:5, letters))
  expect_equal(range_rcpp(y)[1:2], range_rcpp(y))
})

test_that("halting", {
  z <- c(0, 0, -2, 2)
  y <- do_range_dbl(z, -1, 1)
  expect_equal(y[1:2], c(-2, 0))
  expect_equal(y[3], c(3))
  y <- do_range_dbl(z, -3, 1)
  expect_equal(y[4], c(4))

  z <- as.integer(z)
  y <- do_range_int(z, -1L, 1L)
  expect_equal(y[1:2], c(-2, 0))
  expect_equal(y[3], c(3))
  y <- do_range_int(z, -3, 1)
  expect_equal(y[4], c(4))
})

test_that("range on empty", {
  expect_warning(range0 <- range_rcpp(double(0)),
                 regexp = "no non-missing arguments to range_rcpp; returning c(Inf, -Inf).",
                 fixed = TRUE)
  expect_equal(range0, c(Inf, -Inf))
  expect_warning(range_rcpp(integer(0)),
                 regexp = "no non-missing arguments to range_rcpp; returning c(Inf, -Inf).",
                 fixed = TRUE)
  expect_warning(range0i <- range_rcpp(integer(0), integer0_range_is_integer = TRUE),
                 regexp = "no non-missing arguments to range_rcpp; returning c(INT_MAX, -INT_MAX).",
                 fixed = TRUE)
  expect_true(is.integer(range0i))
  expect_true(is.unsorted(range0i))
})

test_that("range on logical", {
  expect_equal(range_rcpp(logical(100)), c(range(logical(100)), 1L, 1L))
  expect_equal(range_rcpp(!logical(100)), c(range(!logical(100)), 1L, 1L))
  expect_equal(range_rcpp(c(FALSE, TRUE, FALSE)), c(0L, 1L, 1L, 2L))
  expect_equal(range_rcpp(c(FALSE, TRUE, NA)), c(0L, 1L, 1L, 2L))
  expect_equal(range_rcpp(c(logical(100), TRUE, NA)), c(0L, 1L, 1L, 101L))
  expect_equal(range_rcpp(c(NA, NA, NA)), as.integer(c(NA, NA, NA, NA)))
})

test_that("range with NA", {
  r1 <- range_rcpp(c(NA, 1L, 2L, 3L))
  expect_equal(r1, c(1, 3, 2, 4))

  r1d <- range_rcpp(as.double(c(NA, 1L, 2L, 3L)))
  expect_equal(r1d, c(1, 3, 2, 4))

  r2 <- range_rcpp(c(NA, NA, 1L, 2L, 3L))
  expect_equal(r2, c(1, 3, 3, 5))

  r2d <- range_rcpp(as.double(c(NA, NA, 1L, 2L, 3L)))
  expect_equal(r2d, c(1, 3, 3, 5))

  r3 <- range_rcpp(c(1L, 9L, -1L, NA, NA))
  expect_equal(r3, c(-1, 9, 3, 2))

  r3d <- range_rcpp(as.double(c(1L, 9L, -1L, NA, NA)))
  expect_equal(r3d, c(-1, 9, 3, 2))

})

test_that("long range", {
  skip_on_cran()
  skip_on_travis()
  x <- tryCatch(integer(.Machine$integer.max + 5),
                error = function(e) {
                  skip("vector specified too large.")
                })
  R <- range_rcpp(x)
  expect_equal(R[1], 0)
  expect_equal(R[2], 0)
})
