test_that("count_logical works", {
  expect_error(count_logical(double(1)), "logical")
  expect_identical(count_logical(c(TRUE, FALSE, NA)), c(1L, 1L, 1L))
  expect_identical(count_logical(c(TRUE, FALSE, TRUE)), c(1L, 2L, 0L))
  is64bit_ <- tryCatch(is.integer(seq_len(2^31 - 1)), error = function(e) FALSE)
  skip_if_not(is64bit_)
  skip_on_travis()
  skip_on_cran()
  skip_if_not_installed("parallel")
  expect_equal(count_logical(logical(1e10)), c(1e10, 0, 0))
  expect_equal(count_logical(logical(1e10), parallel::detectCores()), c(1e10, 0, 0))
})


test_that("count_logical long", {
  skip_on_travis()
  skip_on_appveyor()
  skip_if_not(is64bit())

  x <- logical(.Machine$integer.max + 1)
  cl <- count_logical(x)
  expect_equal(cl, c(.Machine$integer.max + 1, 0, 0))
  skip_on_cran()
  cl <- count_logical(x, nThread = 2)
  expect_equal(cl, c(.Machine$integer.max + 1, 0, 0))
  x <- NULL
  x <- allocate0_except(.Machine$integer.max + 1, 2, NA_integer_, nThread = 2)
  storage.mode(x) <- "logical"
  cl <- count_logical(x, nThread = 2)
  expect_equal(cl, c(.Machine$integer.max, 0, 1))
})
