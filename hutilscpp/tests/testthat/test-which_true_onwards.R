context("test-which_true_onwards")

test_that("Error handling", {
  expect_error(which_true_onwards(1), "logical")
  expect_error(which_true_onwards(NA), "missing")
})

test_that("which_true_onwards works", {
  expect_identical(which_true_onwards(c(TRUE, TRUE, FALSE, TRUE, TRUE)),
                   4L)
  expect_identical(which_true_onwards(c(TRUE, FALSE)),
                   0L)
  expect_identical(which_true_onwards(c(TRUE, TRUE)),
                   1L)
  expect_identical(which_true_onwards(c(FALSE)),
                   0L)
})
