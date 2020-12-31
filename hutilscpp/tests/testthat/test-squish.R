context("test-squish")

test_that("squish works", {
  x <- c(0, 1, 1.5, 2.2)
  expect_equal(max(squish(x, 0, 2)), 2)
  expect_equal(max(squish(x, 0, 3)), max(x))
  expect_equal(min(squish(x, 0, 3)), 0)
  expect_equal(min(squish(x, -1, 3)), 0)
  expect_equal(min(squish(x, 1.5, 3)), 1.5)
  xi <- sample(1:10)
  expect_identical(squish(xi, 1L, 10L), xi)
  expect_identical(squish(xi, 1L, 1L), rep_len(1L, 10L))
  expect_identical(squish(xi, 1L, 1L), rep_len(1L, 10L))
  expect_identical(min(squish(xi, 5L, 6L)), 5L)
})

test_that("Error handling", {
  expect_error(squish(1:5, 1:2, 1),
               regexp = "`length(a) = 2`, but must be length-one.",
               fixed = TRUE)
  expect_error(squish(1:5, 1, integer(0)),
               regexp = "`length(b) = 0`, but must be length-one.",
               fixed = TRUE)
  expect_error(squish(1:5, 1L, 2L, in_place = "foo"),
               regexp = "`in_place` was type character but must be logical.",
               fixed = TRUE)
  expect_error(squish(1:5, 1, 2), regexp = "type integer")
  expect_error(squish("a", "a", "b"), "`x` was type character", fixed = TRUE)
})

test_that("Corner cases", {
  expect_identical(squish(NULL), NULL)
})
