test_that("do_which_in works", {
  x <- 1:10
  y <- 2:5
  expect_equal(do_which_in(x, y), which(x %in% y))
  expect_equal(do_which_in(y, x), which(y %in% x))
  x <- c(NA, x)
  expect_equal(do_which_in(x, y), which(x %in% y))
})

test_that("error handling", {
  skip_on_travis()
  skip_on_cran()
  skip_if_not(is64bit())
  skip_on_ci()
  skip_on_covr()
  expect_error(do_which_in(integer(.Machine$integer.max + 1), 1:10), "exceeds")
})
