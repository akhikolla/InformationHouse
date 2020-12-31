fpath = system.file("testdata",
                    "ploglik_cmodStd_r_psill.rda",
                    package = "gear")
load(fpath)

scmod = cmod_std(model = "exponential", psill = 1, r = 1,
                 angle = geoR_ml$aniso.pars[1],
                 ratio = 1/geoR_ml$aniso.pars[2],
                 radians = TRUE, invert = TRUE)
weights = rep(1, length(y))
gear_ml = optimx::optimx(par = c(1, 1), 
                         fn = ploglik_cmodStd_r_psill, 
                         lower = c(0.001, 0.001),
                         upper = c(50, 50),
                         method = "L-BFGS-B", 
                         x = x, y = y, d = d, nugget = geoR_ml$nugget, 
                         weights = weights, scmod = scmod,
                         control = list(dowarn = FALSE))
gear_reml = optimx::optimx(par = c(1, 1), 
                           fn = ploglik_cmodStd_r_psill, 
                           lower = c(0.001, 0.001),
                           upper = c(50, 50),
                           method = "L-BFGS-B", 
                           x = x, y = y, d = d, nugget = geoR_ml$nugget, 
                           weights = weights, scmod = scmod,
                           reml = TRUE,
                           control = list(dowarn = FALSE))
test_that("ploglik_cmodStd_r_psill accuracy (geoR)", {
  expect_true(max(abs(geoR_ml$cov.pars - c(gear_ml$p2, gear_ml$p1/geoR_ml$aniso.pars[2]))) < 1e-2)
  expect_equal(geoR_ml$loglik, gear_ml$value/-2)
  expect_true(max(abs(geoR_reml$cov.pars - c(gear_reml$p2, gear_reml$p1/geoR_reml$aniso.pars[2]))) < 1e-4)
  expect_equal(geoR_reml$loglik, gear_reml$value/-2)
})

cmod = scmod
cmod$evar = geoR_ml$nugget
data = data.frame(y = y, x1 = coords[,1], x2 = coords[,2])
object = geolm(y ~ 1, data = data, coordnames = c("x1", "x2"),
               mod = cmod)
object_ml = estimate(object, method = "L-BFGS-B",
                     upper = list(psill = 50, r = 50),
                     est_nugget = FALSE, est_par3 = FALSE,
                     est_angle = FALSE, est_ratio = FALSE,
                     verbose = FALSE)
object_reml = estimate(object, method = "L-BFGS-B",
                       upper = list(psill = 50, r = 50),
                       est_nugget = FALSE, est_par3 = FALSE,
                       est_angle = FALSE, est_ratio = FALSE,
                       verbose = FALSE, reml = TRUE)

test_that("estimate r_psill accuracy (geoR)", {
  expect_true(abs(geoR_ml$cov.pars[1] - object_ml$mod$psill) < 1e-2)
  expect_true(abs(geoR_ml$cov.pars[2] - object_ml$mod$r/geoR_ml$aniso.pars[2]) < 1e-2)
  expect_true(abs(geoR_ml$beta - object_ml$coeff) < 1e-2)
  expect_true(abs(geoR_ml$loglik - object_ml$loglik) < 1e-4)
  expect_equivalent(gear_ml[1:3], object_ml$optimx[1:3])
  expect_true(abs(geoR_reml$cov.pars[1] - object_reml$mod$psill) < 1e-2)
  expect_true(abs(geoR_reml$cov.pars[2] - object_reml$mod$r/geoR_reml$aniso.pars[2]) < 1e-2)
  expect_true(abs(geoR_reml$beta - object_reml$coeff) < 1e-2)
  expect_true(abs(geoR_reml_noaniso_loglik - object_reml$loglik) < 1e-1/2)
  expect_equivalent(gear_reml[1:3], object_reml$optimx[1:3])
})
