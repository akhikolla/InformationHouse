test_that("rev works", {
  expect_equal(hutilscpp_rev(integer(0)), integer(0))
  expect_equal(hutilscpp_rev(1:5), 5:1)
  expect_equal(hutilscpp_rev(1:5 + 0), 5:1)
  expect_equal(hutilscpp_rev(c(1, 0, 2)), c(2, 0, 1))
})
