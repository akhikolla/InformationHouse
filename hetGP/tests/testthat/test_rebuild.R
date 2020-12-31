library(hetGP)
context("rebuild")

test_that("rebuild",{
  
  # 1D test
  set.seed(32)
  ## motorcycle data
  library(MASS)
  X <- matrix(mcycle$times, ncol = 1)
  Z <- mcycle$accel
  ## Model fitting
  model <- mleHetGP(X = X, Z = Z, lower = 0.1, upper = 50)

  # Remove internal elements, e.g., to save it
  model1 <- strip(model)
  
  # Rebuild as initial
  model1 <- rebuild(model, robust = T)
  
  xgrid <- matrix(seq(0, 60, length.out = 301), ncol = 1) 

  p0 <- predict(model, xgrid)
  p1 <- predict(model1, xgrid)
  
  expect_equal(p0$mean, p1$mean, tol = 1e-8)
  expect_equal(p0$sd2, p1$sd2, tol = 1e-8)
  expect_equal(p0$nugs, p1$nugs, tol = 1e-8)
  
  # Same for hetTP
  model <- mleHetTP(X = X, Z = Z, lower = 0.1, upper = 50)
  
  # Remove internal elements, e.g., to save it
  model1 <- strip(model)
  
  # Rebuild as initial
  model1 <- rebuild(model, robust = T)
  
  xgrid <- matrix(seq(0, 60, length.out = 301), ncol = 1) 
  
  p0 <- predict(model, xgrid)
  p1 <- predict(model1, xgrid)
  
  expect_equal(p0$mean, p1$mean, tol = 1e-8)
  expect_equal(p0$sd2, p1$sd2, tol = 1e-8)
  expect_equal(p0$nugs, p1$nugs, tol = 1e-8)
  
  
  ##############################################################################
  ## 2D
  set.seed(1)
  nvar <- 2
  
  ## Branin redefined in [0,1]^2
  branin <- function(x){
    if(is.null(nrow(x)))
      x <- matrix(x, nrow = 1)
    x1 <- x[,1] * 15 - 5
    x2 <- x[,2] * 15
    (x2 - 5/(4 * pi^2) * (x1^2) + 5/pi * x1 - 6)^2 + 10 * (1 - 1/(8 * pi)) * cos(x1) + 10
  }
  
  ## Noise field via standard deviation
  noiseFun <- function(x){
    if(is.null(nrow(x)))
      x <- matrix(x, nrow = 1)
    return(1/5*(3*(2 + 2*sin(x[,1]*pi)*cos(x[,2]*3*pi) + 5*rowSums(x^2))))
  }
  
  ## data generating function combining mean and noise fields
  ftest <- function(x){
    return(branin(x) + rnorm(nrow(x), mean = 0, sd = noiseFun(x)))
  }
  
  ## Grid of predictive locations
  ngrid <- 51
  xgrid <- matrix(seq(0, 1, length.out = ngrid), ncol = 1) 
  Xgrid <- as.matrix(expand.grid(xgrid, xgrid))
  
  ## Unique (randomly chosen) design locations
  n <- 50
  Xu <- matrix(runif(n * 2), n)
  X <- Xu[sample(1:n, 20*n, replace = TRUE),]
  
  ## obtain training data response at design locations X
  Z <- ftest(X)
  
  ## Formating of data for model creation (find replicated observations) 
  prdata <- find_reps(X, Z, rescale = FALSE, normalize = FALSE)
  
  ## Model fitting
  model <- mleHetGP(X = list(X0 = prdata$X0, Z0 = prdata$Z0, mult = prdata$mult), Z = prdata$Z,
                    lower = rep(0.01, nvar), upper = rep(10, nvar),
                    covtype = "Matern5_2")
  
  # Remove internal elements, e.g., to save it
  model1 <- strip(model)
  
  # Rebuild as initial
  model1 <- rebuild(model)
  
  p0 <- predict(model, Xgrid)
  p1 <- predict(model1, Xgrid)
  
  expect_equal(p0$mean, p1$mean, tol = 1e-8)
  expect_equal(p0$sd2, p1$sd2, tol = 1e-8)
  expect_equal(p0$nugs, p1$nugs, tol = 1e-8)
  
  # test after update
  Xnew <- matrix(runif(2), 1)
  Znew <- ftest(Xnew)
  
  model <- update(model, Xnew = Xnew, Znew = Znew)
  model1 <- update(model1, Xnew = Xnew, Znew = Znew)
  
  p0 <- predict(model, Xgrid)
  p1 <- predict(model1, Xgrid)
  
  expect_equal(p0$mean, p1$mean, tol = 1e-8)
  expect_equal(p0$sd2, p1$sd2, tol = 1e-8)
  expect_equal(p0$nugs, p1$nugs, tol = 1e-8)
})
