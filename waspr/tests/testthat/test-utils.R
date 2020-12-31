context("utils")

test_that("check mode and hpd function",{

  vec = c(1,2,4,4,4,4,4,4,5,7)

  expect_equal(mode_est(vec), 4)
  expect_equal(hpd_est(vec), c(1,7))

})
