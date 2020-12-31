context("test-cumsum_reset")

test_that("Error handling", {
  expect_error(cumsum_reset(1),
               regexp = "logical")
  expect_error(cumsum_reset(c(NA, TRUE), "missing value"))
  expect_error(cumsum_reset(1:5 > 0, 1:6),
               "same lengths")
  expect_error(cumsum_reset(1:5, "a"))
  expect_error(cumsum_reset(1:26 > 0, letters))
})

test_that("cumsum_reset works", {
  expect_identical(cumsum_reset(c(TRUE, TRUE, FALSE, TRUE, TRUE)),
                   c(1L, 2L, 0L, 1L, 2L))
  expect_identical(cumsum_reset(x = c(TRUE, TRUE, FALSE, TRUE, TRUE),
                                y = 1:5),
                   c(1L, 3L, 0L, 4L, 9L))
  expect_identical(cumsum_reset(x = c(TRUE, TRUE, FALSE, TRUE, TRUE),
                                y = as.double(1:5)),
                   c(1, 3, 0, 4, 9))
  expect_identical(cumsum_reset(FALSE), 0L)
  expect_identical(cumsum_reset(FALSE, 1L), 0L)
  expect_identical(cumsum_reset(FALSE, 2), 0)
})

test_that("Corner cases", {
  expect_identical(cumsum_reset(logical(0)), integer(0L))
  expect_identical(cumsum_reset(logical(0), double(0)), double(0L))
})
