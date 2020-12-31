# save output
fpath = system.file("testdata",  "ll_xycholv_geoR.rda", 
                    package = "gear")
load(fpath)

gear_ll = ll_xycholv(x, y, cholv, reml = FALSE, minus2 = FALSE)
gear_ll_reml = ll_xycholv(x, y, cholv, reml = TRUE, minus2 = FALSE)

test_that("check accuracy of ll_xycholv (geoR)", {
  expect_equal(geoR_ll, gear_ll)
  expect_equal(geoR_ll_reml, gear_ll_reml)
})

# internal calculations to determine mu
cholvix = backsolve(cholv, x, transpose = TRUE)
vix = forwardsolve(cholv, cholvix, upper.tri = TRUE)
xtvix = crossprod(cholvix)
cholxtvix = chol(xtvix)
map2coeff = solve_chol(cholxtvix, t(vix))
mu = coeff = c(map2coeff %*% y)

gear_ll_mu = ll_xycholv(x = NULL, y, cholv, mu = mu,
                        reml = FALSE, minus2 = FALSE)

test_that("check accuracy of ll_xycholv with !is.null(mu)", {
  expect_equal(gear_ll, gear_ll_mu)
})
