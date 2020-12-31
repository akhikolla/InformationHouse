context("test-are_even")

test_that("even works", {
  expect_equal(are_even(c(0, 1, 4)),
               c(TRUE, FALSE, TRUE))
  expect_equal(are_even(c(0L, 1L, 4L)),
               c(TRUE, FALSE, TRUE))
  expect_equal(are_even(-c(0L, 1L, 4L)),
               c(TRUE, FALSE, TRUE))
  expect_equal(which_are_even(c(0, 1, 4)),
               c(1L, 3L))
  expect_identical(which_are_even(c(0L, 1L, 4L)),
                   c(1L, 3L))
  expect_identical(which_are_even(c(5L, 1L, 3L)),
                   integer(0L))
})

test_that("error handling", {
  expect_warning(are_even(c(0, 0.5, 1)),
                 regexp = "`x` was type double, but element 2 = 0.5 was not an integer value",
                 fixed = TRUE)
  expect_warning(which_are_even(c(0, 0.5, 1)),
                 regexp = "`x` was type double, but element 2 = 0.5 was not an integer value",
                 fixed = TRUE)
  expect_error(are_even("a"),
               "`x` was not an integer or double.",
               fixed = TRUE)
  expect_error(which_are_even("a"),
               "`x` was not an integer or double.",
               fixed = TRUE)

})

test_that("long vectors", {
  skip_on_cran()
  skip_on_appveyor()
  skip_on_travis()
  expect_error(which_are_even(integer(1e10)),
               # long for 64 bit, large for 32bit
               regexp = "long|large")
})

test_that("non-finite values", {
  expect_warning(are_even(c(NA, 1)))
  expect_warning(are_even(c(NaN, 1)))
  expect_warning(are_even(c(Inf, 1)))
  expect_warning(are_even(c(-Inf, 1)))

  expect_warning(which_are_even(c(NA, 1)))
  expect_warning(which_are_even(c(NaN, 1)))
  expect_warning(which_are_even(c(Inf, 1)))
  expect_warning(which_are_even(c(-Inf, 1)))

  expect_warning(wh <- which_are_even(c(NA, 1, 2, Inf, -Inf, 7, 8, 8)))
  expect_equal(wh, which((c(NA, 1, 2, Inf, -Inf, 7, 8, 8) %% 2) == 0))
})

test_that("NA", {
  expect_equal(do_are_even(c(NA, 1L, 2L), double(0), 0L), c(NA, FALSE, TRUE))
  expect_equal(do_are_even(integer(0), c(NA, 1, 2), 0L), c(NA, FALSE, TRUE))
})

test_that("large", {
  skip_if_not_installed("withr")
  withr::with_options(list(warn = -1), {
    expect_true(all(are_even(c(0L, 10L, 1e10))))
  })
})

test_that("keep_nas", {
  expect_true(all(are_even(c(0L, 10L), nThread = 1L, keep_nas = FALSE)))
  skip_on_cran()
  expect_false(any(are_even(c(11L, 13L), nThread = 2L, keep_nas = FALSE)))
})
