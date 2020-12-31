data = data.frame(y = rnorm(10), x1 = runif(10),
                 x2 = runif(10))
d = as.matrix(dist(data[,c("x1", "x2")]))
coordnames = ~ x1 + x2
n = nrow(d)
mod = cmod_man(v = exp(-d) + diag(n), evar = 1)
geolm_man = geolm(y ~ x1, data = data, mod = mod,
                  coordnames = coordnames)
mod2 = cmod_std(model = "exponential", psill = 1, r = 1, 
                fvar = 1)
geolm_std = geolm(y ~ x1, data = data, mod = mod2,
                  coordnames = coordnames)

test_that("geolm_fit.cmodMan loglik is accurate (mu = NULL)", {
  expect_equivalent(geolm_man$loglik, geolm_std$loglik)
})

geolm_man2 = geolm(y ~ x1, data = data, mod = mod, mu = 1,
                   coordnames = coordnames)
geolm_std2 = geolm(y ~ x1, data = data, mod = mod2, mu = 1,
                   coordnames = coordnames)

test_that("geolm_fit.cmodMan loglik is accurate (mu)", {
  expect_equivalent(geolm_man2$loglik, geolm_std2$loglik)
})

