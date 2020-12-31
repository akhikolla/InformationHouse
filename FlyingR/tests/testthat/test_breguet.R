context("data input in method breguet")

test_that("Empty fat mass throws an error", {
  expect_error(.breguet(3.77, 1.6, 0, 2, 0.3))
})


# test_that("Factor in ordo other 1 and 2 throws error", {
#   expect_error(.breguet(c(3.77, 3.33), c(1.6, 1.4), c(0.12, 0.11), c(2,3),
#                         c(0.3, 0.2)))
# })
