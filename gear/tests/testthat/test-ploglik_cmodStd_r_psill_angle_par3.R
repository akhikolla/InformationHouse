fpath = system.file("testdata",
                    "ploglik_cmodStd_r_psill_angle_par3.rda",
                    package = "gear")
load(fpath)

scmod = cmod_std(model = "matern", psill = 1, r = 1, par3 = 1,
                 angle = pi/8, 
                 ratio = 1/geoR_ml$aniso.pars[2],
                 invert = TRUE, radians = TRUE)
weights = rep(1, length(y))
gear_ml = optimx::optimx(par = c(1, 1, pi/8, 1), 
                         fn = ploglik_cmodStd_r_psill_angle_par3, 
                         lower = c(0.01, 0.01, 0, 0.001),
                         upper = c(5, 5, pi - 0.001, 2.5),
                         method = "L-BFGS-B", 
                         x = x, y = y, d = d, nugget = geoR_ml$nugget, 
                         weights = weights, scmod = scmod,
                         reml = FALSE,
                         control = list(dowarn = FALSE))
gear_reml = optimx::optimx(par = c(1, 1, pi/8, 1), 
                           fn = ploglik_cmodStd_r_psill_angle_par3, 
                           lower = c(0.001, 0.001, 0, 0.001),
                           upper = c(5, 5, pi - 0.001, 2.5),
                           method = "L-BFGS-B", 
                           x = x, y = y, d = d, nugget = geoR_reml$nugget, 
                           weights = weights, scmod = scmod,
                           reml = TRUE,
                           control = list(dowarn = FALSE))

test_that("ploglik_cmodStd_r_psill_angle_par3 accuracy (geoR)", {
  expect_true(abs(geoR_ml$cov.pars[1] - gear_ml$p2) < 1e-2)
  expect_true(abs(geoR_ml$cov.pars[2] - gear_ml$p1/geoR_ml$aniso.pars[2] ) < 1e-4)
  expect_true(abs(geoR_ml$aniso.pars[1] - gear_ml$p3) < 1e-2)
  expect_true(abs(geoR_ml$kappa - gear_ml$p4) < 1e-2)
  expect_equal(geoR_ml$loglik, gear_ml$value/-2)
  expect_true(abs(geoR_reml$cov.pars[1] - gear_reml$p2) < 1e-2)
  expect_true(abs(geoR_reml$aniso.pars[1] - gear_reml$p3) < 1e-2)
  expect_true(abs(geoR_reml$kappa - gear_reml$p4) < 1e-2)
  expect_equal(geoR_reml$loglik, gear_reml$value/-2)
})

cmod = scmod
cmod$evar = geoR_reml$nugget

data = data.frame(y = y, x1 = coords[,1], x2 = coords[,2])
object = geolm(y ~ 1, data = data, coordnames = c("x1", "x2"),
               mod = cmod)
object_ml = estimate(object, method = "L-BFGS-B",
                     lower = list(r = 0.01, psill = 0.01, angle = 0, ratio = 0.001),
                     upper = list(r = 5, psill = 5, angle = pi - .001, par3 = 2.5),
                     est_nugget = FALSE, est_par3 = TRUE,
                     est_angle = TRUE, est_ratio = FALSE,
                     verbose = FALSE)
object_reml = estimate(object, method = "L-BFGS-B",
                       lower = list(r = 0.001, psill = 0.001, angle = 0, ratio = 0.001),
                       upper = list(r = 5, psill = 5, angle = pi - .001, par3 = 2.5),
                       est_nugget = FALSE, est_par3 = TRUE,
                       est_angle = TRUE, est_ratio = FALSE,
                       verbose = FALSE, reml = TRUE)

test_that("estimate r_psill_angle_ratio_par3 accuracy (geoR)", {
  expect_true(abs(geoR_ml$cov.pars[1] - object_ml$mod$psill) < 1e-2)
  expect_true(abs(geoR_ml$cov.pars[2] - object_ml$mod$r*object_ml$mod$ratio) < 1e-4)
  expect_true(abs(geoR_ml$aniso.pars[1] - object_ml$mod$angle) < 1e-3)
  expect_true(abs(geoR_ml$aniso.pars[2] - 1/object_ml$mod$ratio) < 1e-2)
  expect_true(abs(geoR_ml$kappa - object_ml$mod$par3) < 1e-2)
  expect_true(abs(geoR_ml_beta - object_ml$coeff) < 1e-2)
  expect_equal(gear_ml$value/-2, object_ml$loglik)
  expect_equivalent(gear_ml[1:5], object_ml$optimx[1:5])
  expect_true(abs(geoR_reml$cov.pars[1] - object_reml$mod$psill) < 1e-2)
  expect_true(abs(geoR_reml$cov.pars[2] - object_reml$mod$r*object_reml$mod$ratio) < 1e-4)
  expect_true(abs(geoR_reml$aniso.pars[1] - object_reml$mod$angle) < 1e-4)
  expect_true(abs(geoR_reml$aniso.pars[2] - 1/object_reml$mod$ratio) < 1e-2)
  expect_true(abs(geoR_reml$kappa - object_reml$mod$par3) < 1e-2)
  expect_true(abs(geoR_reml_beta - object_reml$coeff) < 1e-2)
  expect_equivalent(gear_reml[1:5], object_reml$optimx[1:5])
  expect_true(abs(geoR_reml_noaniso_loglik - object_reml$loglik) < 1e-1)
})

cmod_radians = scmod
cmod_radians$evar = geoR_reml$nugget
cmod_radians$angle = scmod$angle * 180/pi
cmod_radians$radians = FALSE
object_radians = geolm(y ~ 1, data = data, coordnames = c("x1", "x2"),
                       mod = cmod_radians)
object_ml_radians = estimate(object_radians, method = "L-BFGS-B",
                             lower = list(r = 0.01, psill = 0.01, angle = 0, par3 = 0.001),
                             upper = list(r = 5, psill = 5, angle = (pi - .001) * 180/pi, par3 = 2.5),
                             est_nugget = FALSE, est_par3 = TRUE,
                             est_angle = TRUE, est_ratio = FALSE,
                             verbose = FALSE)

test_that("estimate w/ and w/o radians r_psill_angle_par3 accuracy", {
  expect_equal(object_ml$mod$angle, object_ml_radians$mod$angle / 180 * pi)
  expect_equal(object_ml$mod$ratio, object_ml_radians$mod$ratio)
  expect_equal(object_ml$mod$psill, object_ml_radians$mod$psill)
  expect_equal(object_ml$mod$r, object_ml_radians$mod$r)
  expect_equal(object_ml$mod$par3, object_ml_radians$mod$par3)
  expect_equal(object_ml$coeff, object_ml_radians$coeff)
  expect_equal(object_ml$loglik, object_ml_radians$loglik)
  object_ml_radians$optimx$xtime = object_ml$optimx$xtime
  expect_equal(object_ml$optimx, object_ml_radians$optimx)
})




