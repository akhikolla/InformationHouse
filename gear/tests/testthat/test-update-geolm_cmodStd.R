set.seed(10)

# generate response
y = rnorm(10)
# generate coordinates
x1 = runif(10); x2 = runif(10)

# data frame for observed data
data = data.frame(y, x1, x2)
coords = cbind(x1, x2)
psill = 2 # partial sill
r = 4 # range parameter
evar = .1 # error variance
fvar = .1 # add finescale variance
# one can't generally distinguish between evar and fvar, but
# this is done for illustration purposes

cmod_std = cmod_std("exponential", psill = psill, r = r, 
                    evar = evar, fvar = fvar)

cmod_std2 = cmod_std("exponential", psill = psill + 1, r = r + .5, 
                    evar = evar + .01, fvar = fvar)

# check geolm update for universal kriging
gear1 = geolm(y ~ x1 + x2, data = data, mod = cmod_std,
              coordnames = c("x1", "x2"))
gear2 = geolm(y ~ x1 + x2, data = data, mod = cmod_std2, 
              coordnames = c("x1", "x2"))
gear2b = update(gear1, cmod_std2)
gear2$call = NULL
gear2b$call = NULL

test_that("update.geolm_cmodStd uk calculations are correct", {
  expect_true(identical(gear2, gear2b))
})

# geolm for simple kriging
gear1 = geolm(y ~ x1 + x2, data = data,
              coordnames = c("x1", "x2"),
              mod = cmod_std, mu = 2)
gear2 = geolm(y ~ x1 + x2, data = data,
              coordnames = c("x1", "x2"),
              mod = cmod_std2, mu = 2)
gear2b = update(gear1, cmod_std2)
gear2$call = NULL
gear2b$call = NULL
test_that("update._cmodgeolmStd sk calculations are correct", {
  expect_true(identical(gear2, gear2b))
})



