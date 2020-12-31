context('C Estimation')
library(activegp)
library(hetGP)

covtypes <- c("Gaussian", "Matern3_2", "Matern5_2")
for (ct in 1:length(covtypes)) {
    test_that(paste(covtypes[ct], " C eigenspace estimation via GP's is consistent"), {
      # Function to discover ASM on
      set.seed(123)
      nvar <- 2
      n <- 200
      r <- 1
      f <- function(x) sin(sum(x))
      true_sub <- c(1,1) / sqrt(2)
      
      # Initial design
      design <- matrix(runif(nvar*n), ncol = nvar)
      response <- apply(design, 1, f)
      model <- mleHomGP(design, response, init = list(theta = rep(1, nvar)),
                        covtype = covtypes[ct]) 
      
      C_hat <- C_GP(model)
      sub_est <- eigen(C_hat$mat)$vectors[,1:r]
      
      expect_true(subspace_dist(sub_est, true_sub) < 5e-3)
    })
}

test_that("C estimation via GP's is consistent", {
  # Function to discover ASM on
  set.seed(123)
  nvar <- 2
  n <- 400
  r <- 1
  f <- function(x) sin(sum(x))
  true_C <- matrix(1/8 * (3 + 2 * cos(2) - cos(4)), nrow = 2, ncol = 2)
  
  # Initial design
  design <- matrix(runif(nvar*n), ncol = nvar)
  response <- apply(design, 1, f)
  model <- mleHomGP(design, response, lower = rep(1e-4, nvar), upper = rep(1,nvar)) 
  
  C_hat <- C_GP(model)
  sub_est <- eigen(C_hat$mat)$vectors[,1:r]
  
  expect_true(norm(C_hat$mat - true_C, 'F') < 1e-4)
})

test_that("C update works", {
  # Function to discover ASM on
  set.seed(123)
  nvar <- 2
  n <- 5
  r <- 1
  f <- function(x) sin(sum(x))
  
  # Initial design
  design <- matrix(runif(nvar*n), ncol = nvar)
  response <- apply(design, 1, f)
  model <- mleHomGP(design, response, lower = rep(1e-4, nvar), upper = rep(0.5,nvar), known = list(g = 1e-4)) 
  
  C_hat <- C_GP(model)
  
  ## First for one new design
  xnew <- matrix(runif(nvar), ncol = nvar)
  ynew <- f(xnew)
  
  C_up <- update(C_hat, xnew, ynew)
  model1 <- update(model, Xnew = xnew, Znew = ynew, maxit = 0)
  C_1 <- C_GP(model1)
  
  expect_true(norm(C_up$mat - C_1$mat) < 1e-8)
  
  C_up2 <- update_C2(C_hat, xnew, ynew)
  expect_true(norm(C_up2 - C_1$mat) < 1e-6)
  
  ## Then for a batch of 3 points
  Xnew <- matrix(runif(nvar*3), ncol = nvar)
  Znew <- apply(Xnew, 1, f)
  
  C_up_b <- update(C_hat, Xnew, Znew)
  modelb <- update(model, Xnew = Xnew, Znew = Znew, maxit = 0)
  C_1_b <- C_GP(modelb)
  
  expect_true(norm(C_up_b$mat - C_1_b$mat) < 1e-8)
  
})

#' C Variance Norm
#' Gives the Frobenius norm of the variance of C , that is of E[(C - E[C])^2] = E[(C - E[C])(C - E[C])]
scvar <- function(Cs) {
  MU <- Reduce(function(i,j)i + j, Cs) / length(Cs)
  DELTA <- lapply(Cs, function(C) C - MU)
  DELTA2 <- lapply(DELTA, function(DELTA) DELTA %*% DELTA)
  SIGMA <- Reduce(function(i,j) i + j, DELTA2) / length(DELTA2)
  return(norm(SIGMA, 'F'))
}

covtypes <- c("Gaussian", "Matern3_2", "Matern5_2")
for (ct in 1:length(covtypes)) {
  test_that(paste(covtypes[ct], " Kernel: Variance of the trace/frob. norm of Cn+1 is correct"),{
    set.seed(1234)
    nvar <- 3
    n <- 20
    f <- function(x)(DiceKriging::hartman3(x) + rnorm(1, sd = 0.1))
    
    # Initial design
    design <- matrix(runif(nvar*n), ncol = nvar)
    response <- apply(design, 1, f)
    model <- mleHomGP(design, response, lower = rep(1e-4, nvar),
                      covtype = covtypes[ct], upper = rep(0.5, nvar), known = list(beta0 = 0)) 
    
    C_hat <- C_GP(model)
    
    samp_trs <- an_trs <- samp_cvs <- an_cvs <- samp_var2 <- an_var2 <- rep(NA, 25)
    for(j in 1:25){ 
      xnew <- matrix(runif(nvar), ncol = nvar)
      
      ynew <- f(xnew)
      ypred <- predict(model, x = xnew)
      
      # sample trace of C(ynew)
      nsamp <- 1e4
      tr_samps <- rep(NA, nsamp)
      
      Cs_samp <- matrix(NA, nsamp, length(C_hat$mat))

      for(i in 1:nsamp){
        ynewi <- rnorm(1, ypred$mean, sqrt(ypred$sd2 + ypred$nugs))
        Ci <- update_C2(C_hat, xnew, ynewi)
        Cs_samp[i,] <- as.vector(Ci)
        tr_samps[i] <- sum(diag(Ci))
      }
      
      Cvsamp <- matrix(diag(var(Cs_samp)),3)
      samp_cvs[j] <- sqrt(sum(Cvsamp^2))
      an_cvs[j] <- sqrt(C_var(C_hat, xnew))
      
      vtr <- C_tr(C_hat, xnew)
      samp_trs[j] <- var(tr_samps)
      an_trs[j] <- vtr

      samp_var2[j] <- scvar(lapply(1:nsamp, function(i) matrix(Cs_samp[i,], nrow = nvar)))
      an_var2[j] <- sqrt(C_var2(C_hat, xnew))
    }
    
    expect_equal(model$nu_hat*an_trs, samp_trs, tolerance = 1e-1)
    expect_equal(model$nu_hat*an_cvs, samp_cvs, tolerance = 1e-1)
    expect_equal(model$nu_hat*an_var2, samp_var2, tolerance = 1e-1)
  })
}

test_that("Estimation for quadratic model works", {
  # Function to discover ASM on is a simple quadratic
  A <- matrix(c(1, -1, 0, -1, 2, -1.5, 0, -1.5, 4), nrow = 3, byrow = TRUE)
  b <- c(1, 4, 9)
  
  # Quadratic function
  ftest <- function(x, sd = 1e-6){
    if(is.null(dim(x))) x <- matrix(x, nrow = 1)
    return(3 + drop(diag(x %*% A %*% t(x)) + x %*% b) + rnorm(nrow(x), sd = sd))
  }
  
  ntrain <- 10000
  design <- 2 * matrix(runif(ntrain * 3), ntrain) - 1
  response <- ftest(design)
  
  C_hat <- C_Q(design, response)
  expect_equal(C_hat, 4/3 * A %*% A + tcrossprod(b), tol = 1e-7)
  
})

