test_that("isnt_integersih", {
  expect_identical(which_isnt_integerish(c(1, 2, 3)), 0L)
  expect_identical(which_isnt_integerish(c(1, 2, 3.1)), 3L)
  expect_identical(which_isnt_integerish(c(1, 2, -3.1)), 3L)
  expect_identical(which_isnt_integerish(c(1, 2, -2.1)), 3L)
  expect_identical(which_isnt_integerish(0L), 0L)
  expect_identical(which_isnt_integerish("a"), 1L)
  expect_identical(which_isnt_integerish(NaN), 1L)
  expect_identical(which_isnt_integerish(c(1, NaN)), 2L)
  expect_identical(which_isnt_integerish(c(2147483648, NaN)), 1L)
  expect_identical(which_isnt_integerish(c(-2147483648, NaN)), 1L)
})

test_that("is_integerish(altrep)", {
  skip_if_not(is64bit())
  skip_on_appveyor()
  x <- -2e9:2e9
  expect_equal(which_isnt_integerish(x), 0L)
  x <- -3e9:2e9
  expect_equal(which_isnt_integerish(x), 1L)
  x <- 2e9:-3e9
  w <- which_isnt_integerish(x)
  expect_gt(w, 0)
  expect_equal(x[w - 1], -.Machine$integer.max)
  x <- -2e9:3e9
  expect_equal(which_isnt_integerish(x), 2e9 + 2 + .Machine$integer.max)
})

