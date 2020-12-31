context("test-utils")

test_that("checkTF works", {
  expect_null(check_TF(TRUE))
  expect_null(check_TF(FALSE))
  expect_error(check_TF(NA), "NA")
  expect_error(check_TF(1:2 > 0), "length")
  expect_error(check_TF(1L), "integer")
})

test_that("isnt_number", {
  expect_false(isnt_number(1))
  expect_false(isnt_number(NA_integer_, na.bad = FALSE))
  expect_true(isnt_number(NA_integer_))
  expect_false(isnt_number(Inf, infinite.bad = FALSE))
  expect_true(isnt_number(Inf))
  x <- integer(0L)
  expect_error(anyOutside(1:5, double(0), 1),
               regexp = "`a` had length 0, but must be length-one",
               fixed = TRUE)
  expect_true(isnt_number(3e9, int.only = TRUE, na.bad = FALSE))
  expect_true(isnt_number(NaN, int.only = TRUE, na.bad = FALSE))
  expect_false(isnt_number(NA_real_, int.only = TRUE, na.bad = FALSE))
  expect_true(isnt_number(-3e9, int.only = TRUE))
})


