test_that("between works", {
  x <- 1:10
  expect_equal(x %(between)% c(1, 10), 1:10 %in% 2:9)
  expect_equal(x %]between[% c(1, 10), 1:10 %in% c(1, 10))
  expect_equal(x %(between)% c(NA, 1L), between(x, NA_integer_, 1L, incbounds = FALSE))
  expect_equal(x %(between)% c(1L, NA), between(x, 1L, NA_integer_, incbounds = FALSE))
  expect_true(all(x %(between)% c(NA, NA)))
  expect_equal(x %]between[% c(NA, 2L), x >= 2L)
  expect_equal(x %]between[% c(2L, NA), x <= 2L)
  expect_true(all(x %]between[% c(NA_integer_, NA)))
})
