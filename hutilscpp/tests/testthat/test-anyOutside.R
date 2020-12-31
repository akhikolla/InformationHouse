context("test-anyoutside")

test_that("Error handling", {
  expect_error(anyOutside(1:10, "a"), "character")
  expect_error(anyOutside(1:10, 3L, "c"), "character")
  expect_error(anyOutside(1:10, 1L, 10L, nas_absent = ""), "logical")
  expect_error(anyOutside(1:10, 1L, 10L, nas_absent = c(TRUE, FALSE)), "length")
  expect_error(anyOutside(1:10, 1, 10),
               regexp = "`x`, `a`, and `b` must be numeric vectors of the same type.",
               fixed = TRUE)
})

test_that("anyOutside works", {
  expect_identical(anyOutside(1:10, 1L, 10L), 0L)
  expect_identical(anyOutside(1:10, 1L, 1L), 2L)
  expect_identical(anyOutside(1:10, 10L, 1L), 1L)
  expect_identical(anyOutside(1:10, 4L, 5L), 1L)
  expect_identical(anyOutside(1:10, 1L, 5L), 6L)
  expect_identical(anyOutside(seq(0, 1, length.out = 5), 0, 1), 0L)
  expect_identical(anyOutside(seq(0, 1, length.out = 5), 0.9, 1), 1L)
  expect_identical(anyOutside(NULL), 0L)

  expect_identical(anyOutside(1:10 + 0, 1, 5), 6L)
  expect_identical(anyOutside(c(1:5, NA_real_), 1, 5), NA_integer_)
  expect_identical(anyOutside(c(1:5, NA_real_), 1, 5, na_is_outside = TRUE), 6L)
  expect_identical(anyOutside(c(1:5, NA_real_), 1, 5, na_is_outside = FALSE), 0L)
})

test_that("NAs", {
  expect_identical(anyOutside(c(NA, 1L), 1L, 2L), NA_integer_)
  expect_identical(anyOutside(c(NA, 1L), 1L, 2L, nas_absent = FALSE), NA_integer_)
  expect_identical(anyOutside(c(NA, 1L), 1L, 2L, nas_absent = FALSE, na_is_outside = FALSE), 0L)
  expect_identical(anyOutside(c(NA, 1L, 4L), 1L, 2L, na_is_outside = FALSE), 3L)
  expect_identical(anyOutside(c(NA, 1L, 4L), 1L, 2L, nas_absent = FALSE, na_is_outside = TRUE), 1L)

  expect_identical(anyOutside(c(1L, 4L, 3L), -1L, 5L, nas_absent = TRUE), 0L)
  expect_identical(anyOutside(c(1L, 4L, 3L), -1L, 5L, nas_absent = TRUE, na_is_outside = FALSE), 0L)

  expect_identical(anyOutside(c(NA, 1L, 9L), 1L, 10L, na_is_outside = NA), NA_integer_)
  expect_identical(anyOutside(c(NA, 1L, 9L), 1L, 10L, na_is_outside = TRUE), 1L)
  expect_identical(anyOutside(c(NA, 1L, 9L), 1L, 10L, na_is_outside = FALSE), 0L)
  expect_identical(anyOutside(c(1:10, NA), 1L, 7L), 8L)

  dbl_1na5 <- as.double(c(1:5, NA, 5))
  expect_identical(anyOutside(dbl_1na5, 1, 5, na_is_outside = FALSE), 0L)
  expect_identical(anyOutside_dbl(dbl_1na5, 1, 5, nas_present = TRUE, na_is_outside = FALSE), 0L)
  expect_identical(anyOutside(dbl_1na5, 1, 5), NA_integer_)
  expect_identical(anyOutside(dbl_1na5, 1, 4), 5L)
  expect_identical(anyOutside(dbl_1na5, 1, 4, nas_absent = FALSE, na_is_outside = FALSE), 5L)
  expect_identical(anyOutside(dbl_1na5, 1, 5, nas_absent = FALSE, na_is_outside = FALSE), 0L)

})
