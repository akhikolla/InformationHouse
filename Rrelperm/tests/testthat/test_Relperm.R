context("test-check")

test_that("error if one input is missing", {

   expect_error(kr2p_ow(0.1,0.1))

})

test_that("error if one input is missing", {

   expect_error(kr2p_gl(0.1,0.1))

})


test_that("error if one input is missing", {

   expect_error(kr3p_Baker(0.1, 0.1))

})


test_that("error if SWCRIT is less than SWCON", {

   expect_error(kr3p_Baker(0.15, 0.1, 0.15, 0.15, 0.2, 0.2, 0.05, 0.05, 0.4, 1, 0.3, 3, 2, 4, 2.5, 101))

})
