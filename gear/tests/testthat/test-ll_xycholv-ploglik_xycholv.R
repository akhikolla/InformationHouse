set.seed(19)
y = rnorm(10)
x = matrix(rep(1, length(y)))
coords = matrix(runif(length(y) * 2), ncol = 2)
d = as.matrix(dist(coords))
pv = exp(-d/3) + 0.1 * diag(length(y))
est_psill = ploglik_xycholv(x, y, chol(pv), return_ll = FALSE)
est_psill_reml = ploglik_xycholv(x, y, chol(pv), reml = TRUE, return_ll = FALSE)
v = pv * est_psill
v_reml = pv * est_psill_reml
# same result
test_that("ploglik_xycholv and ll_xycholv produce same results", {
  expect_equal(ploglik_xycholv(x, y, chol(pv)), 
               ll_xycholv(x, y, chol(v)))
  expect_equal(ploglik_xycholv(x, y, chol(pv), minus2 = FALSE, reml = TRUE), 
               ll_xycholv(x, y, chol(v_reml), minus2 = FALSE, reml = TRUE))
})

# internal calculations to determine mu
cholv = chol(pv)
cholvix = backsolve(cholv, x, transpose = TRUE)
vix = forwardsolve(cholv, cholvix, upper.tri = TRUE)
xtvix = crossprod(cholvix)
cholxtvix = chol(xtvix)
map2coeff = solve_chol(cholxtvix, t(vix))
mu = coeff = c(map2coeff %*% y)

test_that("ploglik_xycholv accuracy with !is.null(mu)", {
  expect_equal(
  ploglik_xycholv(x = NULL, y, cholv, mu = mu),
  ploglik_xycholv(x, y, cholv))
})


                
