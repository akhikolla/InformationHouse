context("test-do_which_first")

test_that("do_which_first works", {
  expect_equal(do_which_first(c(FALSE, FALSE)), 0)
  expect_equal(do_which_first(c(FALSE, TRUE)), 2)
  expect_equal(do_which_first(c(FALSE, TRUE, TRUE)), 2)
  expect_equal(do_which_last(c(FALSE, TRUE, TRUE)), 3)
  expect_equal(do_which_last(c(TRUE, FALSE, FALSE)), 1)
  expect_equal(do_which_last(c(FALSE, FALSE, FALSE)), 0)
  expect_equal(do_which_first_false(c(TRUE, TRUE, FALSE)), 3)
  expect_equal(do_which_first_false(c(FALSE, TRUE, FALSE)), 1)
  expect_equal(do_which_first_false(c(TRUE)), 0)

})
