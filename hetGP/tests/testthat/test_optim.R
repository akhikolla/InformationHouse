library(hetGP)
library(numDeriv)
context("optim")

test_that("optim",{
  
  # Start with gradient predictions
  
  ##------------------------------------------------------------
  ## Example 2: 2D Heteroskedastic GP modeling
  ##------------------------------------------------------------
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
    return(10*(3*(2 + 2*sin(x[,1]*pi)*cos(x[,2]*3*pi) + 5*rowSums(x^2))))
  }
  
  ## data generating function combining mean and noise fields
  ftest <- function(x){
    return(branin(x) + rnorm(nrow(x), mean = 0, sd = noiseFun(x)))
  }
  
  ## Grid of predictive locations
  ngrid <- 31
  xgrid <- matrix(seq(0, 1, length.out = ngrid), ncol = 1) 
  Xgrid <- as.matrix(expand.grid(xgrid, xgrid))
  
  ## Unique (randomly chosen) design locations
  n <- 200 # 50
  Xu <- matrix(runif(n * 2), n)
  
  ## Select replication sites randomly
  X <- Xu#[sample(1:n, 20*n, replace = TRUE),]
  
  ## obtain training data response at design locations X
  Z <- ftest(X)
  
  ## Formating of data for model creation (find replicated observations) 
  prdata <- find_reps(X, Z, rescale = FALSE, normalize = FALSE)
  
  ## Model fitting
  model <- mleHomGP(X = list(X0 = prdata$X0, Z0 = prdata$Z0, mult = prdata$mult), Z = prdata$Z,
                    lower = rep(0.01, nvar), upper = rep(10, nvar),  known = list(beta0 = 0),
                    covtype = "Matern5_2")
  
  predictions <- predict(x = Xgrid, object =  model)
  # for finite difference
  eps <- 1e-6
  predictions1 <- predict(x = Xgrid + matrix(c(eps, 0), nrow(Xgrid), 2, byrow = T), object =  model)
  predictions2 <- predict(x = Xgrid + matrix(c(0, eps), nrow(Xgrid), 2, byrow = T), object =  model)
  
  # pfm <- function(x){predict(model, matrix(x, nrow = 1))$mean}
  # pfs2 <- function(x){predict(model, matrix(x, nrow = 1))$sd2}
  # 
  # grad_refm <- t(apply(Xgrid, 1, grad, func = pfm, method.args = list(d = 6)))
  # grad_refs2 <- t(apply(Xgrid, 1, grad, func = pfs2, method.args = list(d = 6)))
  
  grad_preds <- hetGP:::predict_gr(model, Xgrid)
  expect_equal(grad_preds$mean, 
               cbind((predictions1$mean - predictions$mean)/eps, (predictions2$mean - predictions$mean)/eps), tol = 1e-4)
  
  ## Estimation of grad_sd2 is less precise
  expect_equal(grad_preds$sd2,
               cbind((predictions1$sd2 - predictions$sd2)/eps, (predictions2$sd2 - predictions$sd2)/eps), tol = 1e-4)
  
  ### Test of deriv_crit_EI
  grads_ei_ref <- t(apply(Xgrid, 1, grad, func = crit_EI, model = model, cst = -30))
  grads_ei <- hetGP:::deriv_crit_EI(Xgrid, model, cst = -30)
  expect_equal(grads_ei, grads_ei_ref, tol = 1e-8)

  
  ## Same for TPs
  ## Model fitting
  model <- mleHomTP(X = list(X0 = prdata$X0, Z0 = prdata$Z0, mult = prdata$mult), Z = prdata$Z,
                    lower = rep(0.01, nvar), upper = rep(10, nvar),  known = list(beta0 = 0),
                    covtype = "Matern5_2")
  
  predictions <- predict(x = Xgrid, object =  model)
  # for finite difference
  eps <- 1e-6
  predictions1 <- predict(x = Xgrid + matrix(c(eps, 0), nrow(Xgrid), 2, byrow = T), object =  model)
  predictions2 <- predict(x = Xgrid + matrix(c(0, eps), nrow(Xgrid), 2, byrow = T), object =  model)

  grad_preds <- hetGP:::predict_gr(model, Xgrid)
  expect_equal(grad_preds$mean, 
               cbind((predictions1$mean - predictions$mean)/eps, (predictions2$mean - predictions$mean)/eps), tol = 1e-4)
  
  expect_equal(grad_preds$sd2,
               cbind((predictions1$sd2 - predictions$sd2)/eps, (predictions2$sd2 - predictions$sd2)/eps), tol = 1e-4)
  
  ### Test of deriv_crit_EI
  grads_ei_ref <- t(apply(Xgrid, 1, grad, func = crit_EI, model = model, cst = -30))
  grads_ei <- deriv_crit_EI(Xgrid, model, cst = -30)
  expect_equal(grads_ei, grads_ei_ref, tol = 1e-8)
})