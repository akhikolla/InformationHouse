context("data input in method breguet_adj")

test_that("Empty fat mass throws an error", {
  expect_error(.breguet_adj(3.77, 1.6, 0, 2, 0.3))
})
