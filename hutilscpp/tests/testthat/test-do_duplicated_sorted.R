test_that("do_duplicateds works", {
  x <- rep(c(-.Machine$integer.max, -1, 4, 2),
           c(3, 4, 5, 6))

  expect_equal(do_duplicated_sorted_dbl(x), duplicated(x))

  x <- as.integer(x)
  expect_equal(do_duplicated_sorted_int(x), duplicated(x))
})

