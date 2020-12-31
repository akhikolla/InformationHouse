# save output
fpath = system.file("testdata",  "loglik_cmodMan_dmvtnorm.rda", package = "gear")
load(fpath)

coordnames = ~ x1 + x2
mod = cmod_man(v = v, evar = 1)
geolm_man = geolm(y ~ x1, data = data, mod = mod,
                  coordnames = coordnames)
geolm_man2 = geolm(y ~ x1, data = data, mod = mod, mu = 1,
                   coordnames = coordnames)

test_that("geolm_fit.cmodMan loglik is accurate", {
  expect_equivalent(geolm_man$loglik, ll_dmvnorm)
  expect_equivalent(geolm_man2$loglik, ll_dmvnorm2)
})