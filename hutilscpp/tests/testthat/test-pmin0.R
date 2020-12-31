test_that("pmin0 works", {
  o <- c(-1, 0, 1)
  expect_equal(pmin0(o), c(-1, 0, 0))
  pmin0(o, in_place = TRUE)
  expect_equal(o, c(-1, 0, 0))

  oi <- -1:5
  expect_warning({oi_ans <- pmin0(oi, in_place = TRUE)})
  expect_equal(oi_ans, pmin(-1:5, 0))
  oi <- copy(oi)
  pmin0(oi, in_place = TRUE)
  expect_equal(max(pmin0(oi)), 0L)
  expect_equal(max(oi), 0L)
  expect_equal(min(oi), min(-1:5))
  expect_identical(max(pmin0(oi, TRUE)), 0L)
  expect_equal(pmin0(integer(0)), integer(0))
})

test_that("Error handling", {
  expect_error(pmin0(list()), "list")
  expect_error(pmin0(""), "numeric")

})

test_that("pmin0 sorted", {
  x <- as.double(seq(-1e6, 1e7, length.out = 3e3))
  expect_identical(pmin0(x), do_pmin0_radix_sorted_dbl(x))
  x <- hutilscpp_rev(x)
  expect_identical(pmin0(x), do_pmin0_radix_sorted_dbl(x))
  x <- hutilscpp_rev(x)
  x <- as.integer(x)
  expect_identical(pmin0(x), do_pmin0_radix_sorted_int(x))
  x <- hutilscpp_rev(x)
  expect_identical(pmin0(x), do_pmin0_radix_sorted_int(x))

  x <- -5:6 + 0L
  pmin0(x, sorted = TRUE, in_place = TRUE)
  expect_equal(x, pmin(-5:6, 0))
  x <- as.double(-5:6) + 0
  expect_equal(pmin0(x, sorted = TRUE), pmin(-5:6, 0))
  x <- 5:-6 + 0L
  pmin0(x, sorted = TRUE, in_place = TRUE)
  expect_equal(x, pmin(5:-6, 0))
  x <- 5:-6 + 0
  pmin0(x, sorted = TRUE, in_place = TRUE)
  expect_equal(x, pmin(5:-6, 0))

  x <- -5:6

  expect_equal(pmin0(x, sorted = TRUE), pmin(-5:6, 0))
  x <- as.double(-5:6)
  expect_equal(pmin0(x, sorted = TRUE), pmin(-5:6, 0))
  x <- 5:-6
  expect_equal(pmin0(x, sorted = TRUE), pmin(5:-6, 0))
  x <- 5:-6 + 0

  expect_equal(pmin0(x, sorted = TRUE), pmin(5:-6, 0))
})

test_that("pmin0 sorted but all negative", {
  expect_equal(pmin0(rep(1L, 10), sorted = TRUE), integer(10))
  expect_equal(pmin0(10:1, sorted = TRUE), integer(10))
  expect_equal(pmin0(rep(1, 10), sorted = TRUE), double(10))
  expect_equal(pmin0(10:1 + 0, sorted = TRUE), double(10))
})

test_that("pmin0 altrep", {
  expect_warning(pmin0(1:10, in_place = TRUE), "ALTREP")
  expect_equal(pmin0(1:10), pmin(1:10, 0L))
  expect_equal(pmin0(-1:-10), pmin(-1:-10, 0L))
  expect_equal(pmin0(-1:10), pmin(-1:10, 0L))
  expect_equal(pmin0(1:-10), pmin(1:-10, 0L))
})

