context("test-as_integer_if_safe")

test_that("as_integer_if_safe works", {
  expect_identical(as_integer_if_safe(4), 4L)
  expect_identical(as_integer_if_safe(4:5), 4:5)
  expect_identical(as_integer_if_safe(4:5 + 0), 4:5)
  expect_identical(as_integer_if_safe(4:5 + 0.5), 4:5 + 0.5)
  expect_identical(as_integer_if_safe(c(1, NA, 2)), c(1L, NA, 2L))
  expect_identical(as_integer_if_safe(c("abc", "def")), c("abc", "def"))
  max_int <- .Machine$integer.max
  expect_identical(as_integer_if_safe(c(-max_int, max_int)),
                   as.integer(c(-max_int, max_int)))
  expect_true(is.double(as_integer_if_safe(c(-max_int - 1, max_int))))
  expect_true(is.double(as_integer_if_safe(c(-max_int, max_int + 1))))
  expect_true(is.double(as_integer_if_safe(c(-max_int, NA, max_int + 1))))
  expect_equal(sum(is.na(as_integer_if_safe(c(-max_int, NA, max_int + 1)))), 1)

  expect_false(is.integer(as_integer_if_safe(-max_int + 0.5)))
  expect_false(is.integer(as_integer_if_safe(-max_int - 0.5)))
  expect_false(is.integer(as_integer_if_safe(c(-max_int, max_int, Inf, 0, 1, 2))))
  expect_equal(as_integer_if_safe(c(NA, NaN, 1, 2, 3)),
               as.integer(c(NA, NaN, 1, 2, 3)))

})

test_that("utils", {
  expect_equal(is_safe2int(Inf), 0)
  expect_equal(is_safe2int(NaN), 2)
  expect_equal(is_safe2int(100), 1)
})

test_that("ubsan warning", {
  # NaN not naturally coercible to int
  expect_true(is.integer(as_integer_if_safe(NaN)))
  expect_true(is.integer(as_integer_if_safe(NA_real_)))
  expect_false(is.integer(as_integer_if_safe(Inf)))
  expect_false(is.integer(as_integer_if_safe(-Inf)))
  expect_error(force_as_integer(double(10), na_code = 0L), regexp = "Internal error")
})

test_that("large even numbers", {
  skip_if_not_installed("withr")
  withr::with_options(list(warn = -1), {
    expect_true(all(are_even(c(.Machine$integer.max + 1, -1e10, 1e10))))
    expect_false(any(are_even(c(.Machine$integer.max + 1, -1e10, 1e10) + 1)))
  })
})

test_that("altrep", {
  skip_if_not(is_altrep(1:100))
  skip_if_not(is64bit())
  expect_true(is.double(as_integer_if_safe(-3e9:3e9)))
})


