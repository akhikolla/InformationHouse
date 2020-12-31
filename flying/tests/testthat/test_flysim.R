
context("Parameter input in flysim")

data("birds")

vec_data <- c(5, 4, 7, 8, 9, 10)

test_that("vector input is not allowed", {
  expect_error(flysim(vec_data))
})

# test_that("wrong method? throw error", {
#   expect_error(flysim(data = birds, method = "bregt"))
# })

data2 <- birds[, -3]
test_that("less number of columns", {
  expect_error(flysim(data = data2))
})

# test_that("Consumption", {
#   expect_error(flysim(data = birds, ctrl = list(consume = 10)))
# })

test_that("alpha input", {
  expect_error(flysim(data = birds, settings = list(alpha = 1)))
})
