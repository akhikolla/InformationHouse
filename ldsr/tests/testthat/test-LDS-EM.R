context("EM")

theta0 <- list(A = matrix(0.5),
               B = matrix(0.5, 1, 7),
               C = matrix(0.5),
               D = matrix(0.5, 1, 7),
               Q = matrix(1),
               R = matrix(1),
               mu1 = matrix(1),
               V1 = matrix(1))
obs <- log(ldsr:::P1annual$Qa)
mu <- mean(obs)
y <- t(obs - mu)
N <- 85

uInst <- vInst <- t(ldsr:::P1pc[322:406])
u <- v <- t(NPpc[601:813])

foreach::registerDoSEQ()

test_that("Numerical results are correct for the instrumental period in the first two iterations for P1", {
  smooth1 <- ldsr:::Kalman_smoother(y, uInst, vInst, theta0)
  theta1 <- ldsr:::Mstep(y, uInst, vInst, smooth1)
  smooth2 <- ldsr:::Kalman_smoother(y, uInst, vInst, theta1)
  theta2 <- ldsr:::Mstep(y, uInst, vInst, smooth2)
  expect_equal(smooth1$lik, -11.678657, tolerance = 1e-6)
  expect_equal(smooth1$X[c(1, 85)], c(1.293356, -0.987671), tolerance = 1e-6)
  expect_equal(theta1$A[1],  0.606066, tolerance = 1e-6)
  expect_equal(theta1$C[1], -0.005995, tolerance = 1e-6)
  expect_equal(theta1$Q[1],  3.640236, tolerance = 1e-6)
  expect_equal(smooth2$lik, -0.114224, tolerance = 1e-6)
  expect_equal(theta2$A[1],  0.603945, tolerance = 1e-6)
  expect_equal(theta2$C[1], -0.012004, tolerance = 1e-6)
  expect_equal(theta2$Q[1],  3.644322, tolerance = 1e-6)
})

test_that("Convergence is correct for the instrumental period", {
  fit <- ldsr:::LDS_EM(y, uInst, vInst, theta0, 100, 1e-5)
  expect_equal(length(fit$liks), 68)
  expect_equal(fit$lik, -0.039093, tolerance = 1e-6)
})

test_that("Fixed restart works", {
  fit <- LDS_reconstruction(NPannual, u, v, start.year = 1800, init = make_init(nrow(u), nrow(v), 2))
  expect_is(fit, "list")
})

test_that("Randomized restarts works", {
  fit <- LDS_reconstruction(NPannual, u, v, start.year = 1800, num.restarts = 2)
  expect_is(fit, "list")
})

test_that("Cross validation works with randomized restarts", {
  cv <- cvLDS(NPannual, u, v, start.year = 1800, num.restarts = 2, Z = make_Z(NPannual$Qa, 2))
  expect_is(cv, "list")
})

test_that("u can be NULL", {
  fit <- LDS_reconstruction(NPannual, u = NULL, v, start.year = 1800, num.restarts = 2)
  expect_is(fit, "list")
  cv <- cvLDS(NPannual, u = NULL, v, start.year = 1800, num.restarts = 2, Z = make_Z(NPannual$Qa, 2))
  expect_is(cv, "list")
})

test_that("v can be NULL", {
  fit <- LDS_reconstruction(NPannual, u, v = NULL, start.year = 1800, num.restarts = 2)
  expect_is(fit, "list")
  cv <- cvLDS(NPannual, u, v = NULL, start.year = 1800, num.restarts = 2, Z = make_Z(NPannual$Qa, 2))
  expect_is(cv, "list")
})

test_that("u and v can have different nrows", {
  fit <- LDS_reconstruction(NPannual, u[1:2, ], v, start.year = 1800, num.restarts = 2)
  expect_is(fit, "list")
  cv <- cvLDS(NPannual, u[1:2, ], v, start.year = 1800, num.restarts = 2, Z = make_Z(NPannual$Qa, 2))
  expect_is(cv, "list")
})
