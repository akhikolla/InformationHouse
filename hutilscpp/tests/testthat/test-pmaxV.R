test_that("pmaxV works", {
  expect_equal(pmaxV(1:5, 2:6), 2:6)
  expect_equal(pmaxV(3:7, 2:6), 3:7)
  expect_equal(pmaxV(c(0, 2), c(1, -1)), c(1, 2))
})

test_that("Error handling", {
  expect_error(pmaxV(1:5, 1:6), "same length")
  expect_error(pmaxV(list(1:2), 1), "list")
  expect_error(pmaxV(1, list(1:2)), "list")
  expect_error(pmaxV("1", 1), "numeric")
  expect_error(pmaxV(1, "1"), "numeric")
})

test_that("in_place = TRUE", {
  x <- copy(1:10)
  pmaxV(x, x + 1L, in_place = TRUE)
  expect_equal(x, 2:11)
})

test_that("swap_xy", {
  expect_equal(pminV(1:5, 2:6 + 0),
               pminV(1:5 + 0, 2:6))
})
