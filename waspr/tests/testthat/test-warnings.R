context("warnings")

test_that("wasp function returns correct errors",{

  #if par.names has wrong length
  expect_error(wasp(pois_logistic, par.names = c("test")))

  #if par.names is not a character vector
  expect_error(wasp(pois_logistic, par.names = c(1,2,3)))

  #if mcmc is not numeric
  expect_error(wasp(as.character(pois_logistic)))

  #if mcmc is not 3-dimensional
  expect_error(wasp(array(rep(1, 16), dim = c(2,2,2,2))))


})

test_that("mode and hpd functions return correct errors",{

  #if x is not numeric
  expect_error(mode_est(as.character(seq(1:100))))
  expect_error(hpd_est(as.character(seq(1:100))))

})
