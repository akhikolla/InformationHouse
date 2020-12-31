fpath = system.file("testdata",
                    "ploglik_cmodStd_r_lambda.rda",
                    package = "gear")
load(fpath)

scmod = cmod_std(model = "exponential", psill = 1, r = 1)
weights = rep(1, length(y))
gear_ml = optimx::optimx(c(1, 1), 
                         fn = ploglik_cmodStd_r_lambda, 
                         lower = c(0.001, 0),
                         upper = c(200, 5),
                         method = "L-BFGS-B", 
                         x = x, y = y, d = d,  
                         weights = weights, scmod = scmod,
                         control = list(dowarn = FALSE))
gear_ml_sigmasq = ploglik_cmodStd_r_lambda(c(gear_ml$p1, gear_ml$p2),
                                           x = x, y = y, d = d, weights = weights,
                                           scmod = scmod, return_ll = FALSE)
gear_reml = optimx::optimx(par = c(1, 1), 
                           fn = ploglik_cmodStd_r_lambda, 
                           lower = c(0.001, 0),
                           upper = c(5000, 5),
                           method = "L-BFGS-B", 
                           x = x, y = y, d = d,  
                           weights = weights, scmod = scmod,
                           reml = TRUE,
                           control = list(dowarn = FALSE))
gear_reml_sigmasq = ploglik_cmodStd_r_lambda(c(gear_reml$p1, gear_reml$p2),
                                             x = x, y = y, d = d, weights = weights,
                                             scmod = scmod, reml = TRUE, 
                                             return_ll = FALSE)
test_that("ploglik_cmodStd_r_lambda accuracy (geoR)", {
  expect_true(abs(geoR_ml$cov.pars[2] - gear_ml$p1) < 1e-5)
  expect_true(abs(geoR_ml$sigmasq - gear_ml_sigmasq) < 1e-5)
  expect_true(abs(geoR_ml$nugget - gear_ml_sigmasq * gear_ml$p2) < 1e-5)
  expect_equal(geoR_ml$loglik, gear_ml$value/-2)
  expect_true(abs(geoR_reml$cov.pars[2] - gear_reml$p1) < 1e-5)
  expect_true(abs(geoR_reml$sigmasq - gear_reml_sigmasq) < 1e-5)
  expect_true(abs(geoR_reml$nugget - gear_reml_sigmasq * gear_reml$p2) < 1e-5)
  expect_equal(geoR_ml$loglik, gear_ml$value/-2)
})

cmod = scmod
cmod$evar = scmod$psill
data = data.frame(y = y, x1 = coords[,1], x2 = coords[,2])
object = geolm(y ~ 1, data = data, coordnames = c("x1", "x2"),
               mod = cmod)
object_ml = estimate(object, method = "L-BFGS-B",
                     lower = list(r = 0.001, lambda = 0),
                     upper = list(r = 200, lambda = 5),
                     est_nugget = TRUE, est_par3 = FALSE,
                     est_angle = FALSE, est_ratio = FALSE,
                     verbose = FALSE)
object_reml = estimate(object, method = "L-BFGS-B",
                       lower = list(r = 0.001, lambda = 0),
                       upper = list(r = 5000, lambda = 5),
                       est_nugget = TRUE, est_par3 = FALSE,
                       est_angle = FALSE, est_ratio = FALSE,
                       verbose = FALSE, reml = TRUE)

test_that("estimate r_lambda_ratio accuracy (geoR)", {
  expect_true(abs(geoR_ml$sigmasq - object_ml$mod$psill) < 1e-1)
  expect_true(abs(geoR_ml$phi - object_ml$mod$r) < 1e-2)
  expect_true(abs(geoR_ml$kappa - object_ml$mod$par3) < 1e-1)
  expect_true(abs(geoR_ml$beta - object_ml$coeff) < 1e-2)
  expect_equal(gear_ml$value/-2, object_ml$loglik)
  expect_equivalent(gear_ml[1:3], object_ml$optimx[1:3])
  expect_true(abs(geoR_reml$sigmasq - object_reml$mod$psill) < 1e-1)
  expect_true(abs(geoR_reml$phi - object_reml$mod$r) < 1e-4)
  expect_true(abs(geoR_reml$kappa - object_reml$mod$par3) < 1e-2)
  expect_true(abs(geoR_reml$beta - object_reml$coeff) < 1e-2)
  expect_equivalent(gear_reml[1:3], object_reml$optimx[1:3])
})
