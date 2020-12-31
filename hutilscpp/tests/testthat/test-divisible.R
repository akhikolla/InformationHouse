test_that("divisible works", {
  expect_equal(divisible(1:100, 3L), (1:100) %% 3L == 0L)
  expect_equal(divisible16(1:100), (1:100) %% 16L == 0L)
})
