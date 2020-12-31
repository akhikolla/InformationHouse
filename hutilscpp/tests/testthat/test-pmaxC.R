test_that("pmaxC works", {
  x <- 1:10
  expect_equal(pmaxC(x, 3L), pmax(x, 3L))
  expect_equal(pmaxC(x, 3), pmax(x, 3))
})

test_that("pmaxC in-place", {
  o <- c(-1, 0, 1)
  e <- c(0.5, 0.5, 1)
  expect_equal(pmaxC(o, 0.5), c(0.5, 0.5, 1))
  expect_equal(o, c(-1, 0, 1))
  pmaxC(o, 0.5, in_place = TRUE)
  expect_equal(o, c(0.5, 0.5, 1))


  oi <- -1:5
  pmax2 <- pmax(-1:5, 2L)
  expect_equal(pmaxC(oi, 2L), pmax(oi, 2L))
  pmaxC(oi, 2L, in_place = TRUE)
  expect_equal(oi, pmax2)
})

test_that("pmaxC error handling", {
  expect_error(pmaxC("", ""), regexp = "numeric")
  expect_error(pmaxC(list(1), 1L), regexp = "atomic")
  expect_error(pmaxC(1:5, 1:2), "length.*2")
  expect_message(pmaxC(1:5, 0.5), "Output is double.")
  expect_error(pmaxC(1:5, 2.5, dbl_ok = FALSE))
})

test_that("pmaxC_real_real", {
  x <- runif(100)
  y <- x[2]
  expect_equal(pminC(x, y), pmin(x, y))
  x <- c(NaN, x)
  expect_equal(pminC(x, y, keep_nas = TRUE), pmin(x, y))
})

test_that("pmaxC_int_real", {
  x <- 1:5
  o <- pmax(x, 2)
  pmaxC(x, 2, in_place = TRUE)
  expect_equal(x, o)
  xd <- 1:5 + 0
  pmaxC(xd, 2L, in_place = TRUE)
  expect_equal(xd, o)



})

