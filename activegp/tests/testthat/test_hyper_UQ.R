context("hyper_UQ")

test_that("homGP Hessian is correct", {
  library(numDeriv)
  library(hetGP)
  
  set.seed(123)
  nvar <- 2
  n <- 30
  r <- 1
  f <- function(x) sin(sum(x))
  true_C <- matrix(1/8 * (3 + 2 * cos(2) - cos(4)), nrow = 2, ncol = 2)
  
  # Initial design
  design <- matrix(runif(nvar*n), ncol = nvar)
  response <- apply(design, 1, f)
  model <- mleHomGP(design, response)
  
  ll <- drop(logLikH(X0 = model$X0, Z0 = model$Z0,Z = model$Z, mult = model$mult, theta = model$theta,
               g = model$g, beta0 = model$beta0))
  
  expect_equal(ll, model$ll)
  
  gr <- grad(logLikH, model$theta, g = model$g, X0 = model$X0, Z0 = model$Z0,
             Z = model$Z, mult = model$mult, beta0 = model$beta0,
             Delta = model$Delta, k_theta_g = model$k_theta_g, theta_g = model$theta_g,
             logN = model$logN, covtype = model$covtype)
  
  gr_ref <- hetGP:::dlogLikHom(X0 = model$X0, theta = model$theta, g = model$g, Z0 = model$Z0,
                               Z = model$Z, mult = model$mult, beta0 = model$beta0, covtype = model$covtype)
  expect_equal(gr, gr_ref[1:nvar], tol = 1e-4)
  # If gradient is ok, hessian should be
  
})

test_that("hetGP Hessian is correct", {
  library(numDeriv)
  library(hetGP)
  
  set.seed(123)
  nvar <- 2
  n <- 30
  r <- 1
  f <- function(x) sin(sum(x))
  true_C <- matrix(1/8 * (3 + 2 * cos(2) - cos(4)), nrow = 2, ncol = 2)
  
  # Initial design
  design <- matrix(runif(nvar*n), ncol = nvar)
  response <- apply(design, 1, f) + rnorm(n, sd = 1e-3)
  model <- mleHetGP(design, response, settings = list(checkHom = FALSE, linkThetas = "joint"))
  
  ll <- drop(logLikH(X0 = model$X0, Z0 = model$Z0,Z = model$Z, mult = model$mult, theta = model$theta,
                     g = model$g, beta0 = model$beta0, Delta = model$Delta, k_theta_g = model$k_theta_g,
                     theta_g = model$theta_g, logN = model$logN))
  
  expect_equal(ll, drop(model$ll_non_pen))
  
  gr <- grad(logLikH, model$theta, g = model$g, X0 = model$X0, Z0 = model$Z0,
             Z = model$Z, mult = model$mult, beta0 = model$beta0,
             Delta = model$Delta, k_theta_g = model$k_theta_g, theta_g = model$theta_g,
             logN = model$logN, covtype = model$covtype)
  
  gr_ref <- hetGP:::dlogLikHet(X0 = model$X0, theta = model$theta, g = model$g, Z0 = model$Z0,
                               Z = model$Z, mult = model$mult, beta0 = model$beta0, covtype = model$covtype,
                               Delta = model$Delta, k_theta_g = model$k_theta_g, theta_g = model$theta_g, logN = model$logN,
                               penalty = FALSE)
  expect_equal(gr, gr_ref[1:nvar], tol = 1e-4)
  # If gradient is ok, hessian should be
  
})

