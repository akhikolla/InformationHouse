context("testcorr")

test_that("Correlation CorrMatCauchy works", {
  x1 <- runif(5)
  x2 <- runif(4)
  th <- c(.05,.9,-.3)
  th <- runif(3,-1,1)

  # First check return_numpara is right
  expect_equal(CGGP_internal_CorrMatCauchy(return_numpara=TRUE), 3)
  
  # Get error when you give in theta of wrong size
  expect_error(CGGP_internal_CorrMatCauchy(x1=x1, x2=x2, theta = c(.1,.1)))

  # Now check correlation
  cauchy1 <- CGGP_internal_CorrMatCauchy(x1=x1, x2=x2, theta=th)
  expect_is(cauchy1, "matrix")
  expect_equal(dim(cauchy1), c(5,4))

  # This just copies the function as written, so not really an independent check.
  #  Should pass by definition. But will help in case we change the correlation
  #  function later, or convert it to Rcpp.
  cauchyfunc <- function(a,b,theta) {
    expLS <- exp(3*theta[1])
    expHE <- exp(3*theta[2])
    alpha = 2*exp(3*theta[3]+2)/(1+exp(3*theta[3]+2))
    diffmat <- outer(a, b, Vectorize(function(aa,bb) abs(aa-bb)))
    h = diffmat/expLS
    halpha = h^alpha
    pow = -expHE/alpha
    (1+halpha)^pow
  }
  cauchy2 <- cauchyfunc(x1, x2, theta=th)
  expect_equal(cauchy1, cauchy2, tol=1e-5)

  # Now check that dC is actually grad of C
  cauchy_C_dC <- CGGP_internal_CorrMatCauchy(x1=x1, x2=x2, theta=th, return_dCdtheta=TRUE)
  eps <- 1e-6
  for (i in 1:3) {
    thd <- c(0,0,0)
    thd[i] <- eps
    numdC <- (CGGP_internal_CorrMatCauchy(x1=x1, x2=x2, theta=th+thd) -
              CGGP_internal_CorrMatCauchy(x1=x1, x2=x2, theta=th-thd)) / (2*eps)
    # numdC <- (-CGGP_internal_CorrMatCauchy(x1=x1, x2=x2, theta=th+2*thd) +
    #           8*  CGGP_internal_CorrMatCauchy(x1=x1, x2=x2, theta=th+thd) +
    #           -8*  CGGP_internal_CorrMatCauchy(x1=x1, x2=x2, theta=th-thd) +
    #             CGGP_internal_CorrMatCauchy(x1=x1, x2=x2, theta=th-2*thd)) / (12*eps)
    # Should be more accurate but was exactly the same
    expect_equal(numdC, cauchy_C_dC$dCdtheta[,(1+4*i-4):(4*i)], info = paste("theta dimension with error is",i))
  }
})

test_that("Correlation CorrMatCauchySQT works", {
  x1 <- runif(5)
  x2 <- runif(4)
  th <- c(.05,.9,-.3)
  th <- runif(3,-1,1)
  
  # First check return_numpara is right
  expect_equal(CGGP_internal_CorrMatCauchySQT(return_numpara=TRUE), 3)
  
  # Get error when you give in theta of wrong size
  expect_error(CGGP_internal_CorrMatCauchySQT(x1=x1, x2=x2, theta = c(.1,.1)))
  
  # Now check correlation
  cauchy1 <- CGGP_internal_CorrMatCauchySQT(x1=x1, x2=x2, theta=th)
  expect_is(cauchy1, "matrix")
  expect_equal(dim(cauchy1), c(5,4))
  
  # This just copies the function as written, so not really an independent check.
  #  Should pass by definition. But will help in case we change the correlation
  #  function later, or convert it to Rcpp.
  cauchyfunc <- function(x1,x2,theta) {
    
    expTILT = exp((theta[3]))
    expLS = exp(3*(theta[1]))
    x1t = (x1+10^(-2))^expTILT
    x2t = (x2+10^(-2))^expTILT
    x1ts = x1t/expLS
    x2ts = x2t/expLS
    
    diffmat =abs(outer(x1ts,x2ts,'-')); 
    expHE = exp(3*(theta[2]))
    h = diffmat
    alpha = 2*exp(5)/(1+exp(5))
    halpha = h^alpha
    pow = -expHE/alpha
    
    (1+halpha)^pow
  }
  cauchy2 <- cauchyfunc(x1, x2, theta=th)
  expect_equal(cauchy1, cauchy2)
  
  # Now check that dC is actually grad of C
  cauchy_C_dC <- CGGP_internal_CorrMatCauchySQT(x1=x1, x2=x2, theta=th, return_dCdtheta=TRUE)
  eps <- 1e-6
  for (i in 1:3) {
    thd <- c(0,0,0)
    thd[i] <- eps
    numdC <- (CGGP_internal_CorrMatCauchySQT(x1=x1, x2=x2, theta=th+thd) -
                CGGP_internal_CorrMatCauchySQT(x1=x1, x2=x2, theta=th-thd)) / (2*eps)
    # Should be more accurate but was exactly the same
    expect_equal(numdC, cauchy_C_dC$dCdtheta[,(1+4*i-4):(4*i)], info = paste("theta dimension with error is",i))
  }
})


test_that("Correlation CorrMatCauchySQ works", {
  x1 <- runif(5)
  x2 <- runif(4)
  # th <- c(.05,.9)
  th <- runif(2,-1,1)
  
  # First check return_numpara is right
  expect_equal(CGGP_internal_CorrMatCauchySQ(return_numpara=TRUE), 2)
  
  # Get error when you give in theta of wrong size
  expect_error(CGGP_internal_CorrMatCauchySQ(x1=x1, x2=x2, theta = c(.1,.1,.4)))
  
  # Now check correlation
  cauchy1 <- CGGP_internal_CorrMatCauchySQ(x1=x1, x2=x2, theta=th)
  expect_is(cauchy1, "matrix")
  expect_equal(dim(cauchy1), c(5,4))
  
  # This just copies the function as written, so not really an independent check.
  #  Should pass by definition. But will help in case we change the correlation
  #  function later, or convert it to Rcpp.
  cauchyfunc <- function(x1,x2,theta) {
    diffmat =abs(outer(x1,x2,'-')); 
    
    expLS = exp(3*theta[1])
    expHE = exp(3*theta[2])
    h = diffmat/expLS
    alpha = 2*exp(0+6)/(1+exp(0+6))
    halpha = h^alpha
    pow = -expHE/alpha
    
    (1-10^(-10))*(1+halpha)^pow+10^(-10)*(diffmat<10^(-4))
  }
  cauchy2 <- cauchyfunc(x1, x2, theta=th)
  expect_equal(cauchy1, cauchy2)
  
  # Now check that dC is actually grad of C
  cauchy_C_dC <- CGGP_internal_CorrMatCauchySQ(x1=x1, x2=x2, theta=th, return_dCdtheta=TRUE)
  eps <- 1e-6
  for (i in 1:2) {
    thd <- c(0,0)
    thd[i] <- eps
    numdC <- (CGGP_internal_CorrMatCauchySQ(x1=x1, x2=x2, theta=th+thd) -
                CGGP_internal_CorrMatCauchySQ(x1=x1, x2=x2, theta=th-thd)) / (2*eps)
    # Should be more accurate but was exactly the same
    expect_equal(numdC, cauchy_C_dC$dCdtheta[,(1+4*i-4):(4*i)], info = paste("theta dimension with error is",i))
  }
})


test_that("Correlation CorrMatGaussian works", {
  x1 <- runif(5)
  x2 <- runif(4)
  th <- runif(1,-1,1)
  
  # First check return_numpara is right
  expect_equal(CGGP_internal_CorrMatGaussian(return_numpara=TRUE), 1)
  
  # Get error when you give in theta of wrong size
  expect_error(CGGP_internal_CorrMatGaussian(x1=x1, x2=x2, theta = c(.1,.1)))
  
  # Now check correlation
  corr1 <- CGGP_internal_CorrMatGaussian(x1=x1, x2=x2, theta=th)
  expect_is(corr1, "matrix")
  expect_equal(dim(corr1), c(5,4))
  
  # This just copies the function as written, so not really an independent check.
  #  Should pass by definition. But will help in case we change the correlation
  #  function later, or convert it to Rcpp.
  gaussianfunc <- function(x1,x2,theta) {
    diffmat =abs(outer(x1,x2,'-')); 
    diffmat2 <- diffmat^2
    expLS = exp(3*theta[1])
    h = diffmat2/expLS
    C = (1-10^(-10))*exp(-h) + 10^(-10)*(diffmat<10^(-4))
    C
  }
  corr2 <- gaussianfunc(x1, x2, theta=th)
  expect_equal(corr1, corr2)
  
  # Now check that dC is actually grad of C
  corr_C_dC <- CGGP_internal_CorrMatGaussian(x1=x1, x2=x2, theta=th, return_dCdtheta=TRUE)
  eps <- 1e-6
  for (i in 1:1) {
    thd <- c(0)
    thd[i] <- eps
    numdC <- (CGGP_internal_CorrMatGaussian(x1=x1, x2=x2, theta=th+thd) -
                CGGP_internal_CorrMatGaussian(x1=x1, x2=x2, theta=th-thd)) / (2*eps)
    # Should be more accurate but was exactly the same
    expect_equal(numdC, corr_C_dC$dCdtheta[,(1+4*i-4):(4*i)], info = paste("theta dimension with error is",i))
  }
})


test_that("Correlation CorrMatMatern32 works", {
  x1 <- runif(5)
  x2 <- runif(4)
  th <- runif(1,-1,1)
  
  # First check return_numpara is right
  expect_equal(CGGP_internal_CorrMatMatern32(return_numpara=TRUE), 1)
  
  # Get error when you give in theta of wrong size
  expect_error(CGGP_internal_CorrMatMatern32(x1=x1, x2=x2, theta = c(.1,.1)))
  
  # Now check correlation
  corr1 <- CGGP_internal_CorrMatMatern32(x1=x1, x2=x2, theta=th)
  expect_is(corr1, "matrix")
  expect_equal(dim(corr1), c(5,4))
  
  # This just copies the function as written, so not really an independent check.
  #  Should pass by definition. But will help in case we change the correlation
  #  function later, or convert it to Rcpp.
  matern32func <- function(x1,x2,theta) {
    diffmat =abs(outer(x1,x2,'-'))
    expLS = exp(3*theta[1])
    h = diffmat/expLS
    C = (1-10^(-10))*(1+sqrt(3)*h)*exp(-sqrt(3)*h) + 10^(-10)*(diffmat<10^(-4))
    C
  }
  corr2 <- matern32func(x1, x2, theta=th)
  expect_equal(corr1, corr2)
  
  # Now check that dC is actually grad of C
  corr_C_dC <- CGGP_internal_CorrMatMatern32(x1=x1, x2=x2, theta=th, return_dCdtheta=TRUE)
  eps <- 1e-6
  for (i in 1:1) {
    thd <- c(0)
    thd[i] <- eps
    numdC <- (CGGP_internal_CorrMatMatern32(x1=x1, x2=x2, theta=th+thd) -
                CGGP_internal_CorrMatMatern32(x1=x1, x2=x2, theta=th-thd)) / (2*eps)
    # Should be more accurate but was exactly the same
    expect_equal(numdC, corr_C_dC$dCdtheta[,(1+4*i-4):(4*i)], info = paste("theta dimension with error is",i))
  }
})



test_that("Correlation CorrMatMatern52 works", {
  x1 <- runif(5)
  x2 <- runif(4)
  th <- runif(1,-1,1)
  
  # First check return_numpara is right
  expect_equal(CGGP_internal_CorrMatMatern52(return_numpara=TRUE), 1)
  
  # Get error when you give in theta of wrong size
  expect_error(CGGP_internal_CorrMatMatern52(x1=x1, x2=x2, theta = c(.1,.1)))
  
  # Now check correlation
  corr1 <- CGGP_internal_CorrMatMatern52(x1=x1, x2=x2, theta=th)
  expect_is(corr1, "matrix")
  expect_equal(dim(corr1), c(5,4))
  
  # This just copies the function as written, so not really an independent check.
  #  Should pass by definition. But will help in case we change the correlation
  #  function later, or convert it to Rcpp.
  matern52func <- function(x1,x2,theta) {
    diffmat =abs(outer(x1,x2,'-'))
    expLS = exp(3*theta[1])
    h = diffmat/expLS
    C = (1-10^(-10))*(1+sqrt(5)*h+5/3*h^2)*exp(-sqrt(5)*h) + 10^(-10)*(diffmat<10^(-4))
    C
  }
  corr2 <- matern52func(x1, x2, theta=th)
  expect_equal(corr1, corr2)
  
  # Now check that dC is actually grad of C
  corr_C_dC <- CGGP_internal_CorrMatMatern52(x1=x1, x2=x2, theta=th, return_dCdtheta=TRUE)
  eps <- 1e-6
  for (i in 1:1) {
    thd <- c(0)
    thd[i] <- eps
    numdC <- (CGGP_internal_CorrMatMatern52(x1=x1, x2=x2, theta=th+thd) -
                CGGP_internal_CorrMatMatern52(x1=x1, x2=x2, theta=th-thd)) / (2*eps)
    # Should be more accurate but was exactly the same
    expect_equal(numdC, corr_C_dC$dCdtheta[,(1+4*i-4):(4*i)], info = paste("theta dimension with error is",i))
  }
})


test_that("Correlation CorrMatPowerExp works", {
  nparam <- 2
  x1 <- runif(5)
  x2 <- runif(4)
  th <- runif(nparam,-1,1)
  
  # First check return_numpara is right
  expect_equal(CGGP_internal_CorrMatPowerExp(return_numpara=TRUE), nparam)
  
  # Get error when you give in theta of wrong size
  expect_error(CGGP_internal_CorrMatPowerExp(x1=x1, x2=x2, theta = rep(0, nparam-1)))
  expect_error(CGGP_internal_CorrMatPowerExp(x1=x1, x2=x2, theta = rep(0, nparam+1)))
  
  # Now check correlation
  corr1 <- CGGP_internal_CorrMatPowerExp(x1=x1, x2=x2, theta=th)
  expect_is(corr1, "matrix")
  expect_equal(dim(corr1), c(5,4))
  
  # This just copies the function as written, so not really an independent check.
  #  Should pass by definition. But will help in case we change the correlation
  #  function later, or convert it to Rcpp.
  PowerExpfunc <- function(x1,x2,theta) {
    diffmat =abs(outer(x1,x2,'-'))
    expLS = exp(3*theta[1])
    minpower <- 1
    maxpower <- 1.95
    alpha <- minpower + (theta[2]+1)/2 * (maxpower - minpower)
    h = diffmat/expLS
    C = (1-10^(-10))*exp(-(h)^alpha) + 10^(-10)*(diffmat<10^(-4))
    C
  }
  corr2 <- PowerExpfunc(x1, x2, theta=th)
  expect_equal(corr1, corr2)
  
  # Now check that dC is actually grad of C
  corr_C_dC <- CGGP_internal_CorrMatPowerExp(x1=x1, x2=x2, theta=th, return_dCdtheta=TRUE)
  eps <- 1e-6
  for (i in 1:nparam) {
    thd <- rep(0, nparam)
    thd[i] <- eps
    numdC <- (CGGP_internal_CorrMatPowerExp(x1=x1, x2=x2, theta=th+thd) -
                CGGP_internal_CorrMatPowerExp(x1=x1, x2=x2, theta=th-thd)) / (2*eps)
    # Should be more accurate but was exactly the same
    expect_equal(numdC, corr_C_dC$dCdtheta[,(1+4*i-4):(4*i)], info = paste("theta dimension with error is",i))
  }
})



test_that("Correlation CorrMatWendland0 works", {
  x1 <- runif(5)
  x2 <- runif(4)
  th <- runif(1,-1,1)
  
  # First check return_numpara is right
  expect_equal(CGGP_internal_CorrMatWendland0(return_numpara=TRUE), 1)
  
  # Get error when you give in theta of wrong size
  expect_error(CGGP_internal_CorrMatWendland0(x1=x1, x2=x2, theta = c(.1,.1)))
  
  # Now check correlation
  corr1 <- CGGP_internal_CorrMatWendland0(x1=x1, x2=x2, theta=th)
  expect_is(corr1, "matrix")
  expect_equal(dim(corr1), c(5,4))
  
  # This just copies the function as written, so not really an independent check.
  #  Should pass by definition. But will help in case we change the correlation
  #  function later, or convert it to Rcpp.
  wendland0func <- function(x1,x2,theta) {
    diffmat =abs(outer(x1,x2,'-'))
    expLS = exp(3*theta[1])
    h = diffmat/expLS
    C = pmax(1 - h, 0)
    C
  }
  corr2 <- wendland0func(x1, x2, theta=th)
  expect_equal(corr1, corr2)
  
  # Now check that dC is actually grad of C
  corr_C_dC <- CGGP_internal_CorrMatWendland0(x1=x1, x2=x2, theta=th, return_dCdtheta=TRUE)
  eps <- 1e-6
  for (i in 1:1) {
    thd <- c(0)
    thd[i] <- eps
    numdC <- (CGGP_internal_CorrMatWendland0(x1=x1, x2=x2, theta=th+thd) -
                CGGP_internal_CorrMatWendland0(x1=x1, x2=x2, theta=th-thd)) / (2*eps)
    # Should be more accurate but was exactly the same
    expect_equal(numdC, corr_C_dC$dCdtheta[,(1+4*i-4):(4*i)], info = paste("theta dimension with error is",i))
  }
})



test_that("Correlation CorrMatWendland1 works", {
  x1 <- runif(5)
  x2 <- runif(4)
  th <- runif(1,-1,1)
  
  # First check return_numpara is right
  expect_equal(CGGP_internal_CorrMatWendland1(return_numpara=TRUE), 1)
  
  # Get error when you give in theta of wrong size
  expect_error(CGGP_internal_CorrMatWendland1(x1=x1, x2=x2, theta = c(.1,.1)))
  
  # Now check correlation
  corr1 <- CGGP_internal_CorrMatWendland1(x1=x1, x2=x2, theta=th)
  expect_is(corr1, "matrix")
  expect_equal(dim(corr1), c(5,4))
  
  # This just copies the function as written, so not really an independent check.
  #  Should pass by definition. But will help in case we change the correlation
  #  function later, or convert it to Rcpp.
  wendland1func <- function(x1,x2,theta) {
    diffmat =abs(outer(x1,x2,'-'))
    expLS = exp(3*theta[1])
    h = diffmat/expLS
    C = pmax(1 - h, 0)^3 * (3*h+1)
    C
  }
  corr2 <- wendland1func(x1, x2, theta=th)
  expect_equal(corr1, corr2)
  
  # Now check that dC is actually grad of C
  corr_C_dC <- CGGP_internal_CorrMatWendland1(x1=x1, x2=x2, theta=th, return_dCdtheta=TRUE)
  eps <- 1e-6
  for (i in 1:1) {
    thd <- c(0)
    thd[i] <- eps
    numdC <- (CGGP_internal_CorrMatWendland1(x1=x1, x2=x2, theta=th+thd) -
                CGGP_internal_CorrMatWendland1(x1=x1, x2=x2, theta=th-thd)) / (2*eps)
    # Should be more accurate but was exactly the same
    expect_equal(numdC, corr_C_dC$dCdtheta[,(1+4*i-4):(4*i)], info = paste("theta dimension with error is",i))
  }
})


test_that("Correlation CorrMatWendland2 works", {
  x1 <- runif(5)
  x2 <- runif(4)
  th <- runif(1,-1,1)
  
  # First check return_numpara is right
  expect_equal(CGGP_internal_CorrMatWendland2(return_numpara=TRUE), 1)
  
  # Get error when you give in theta of wrong size
  expect_error(CGGP_internal_CorrMatWendland2(x1=x1, x2=x2, theta = c(.1,.1)))
  
  # Now check correlation
  corr1 <- CGGP_internal_CorrMatWendland2(x1=x1, x2=x2, theta=th)
  expect_is(corr1, "matrix")
  expect_equal(dim(corr1), c(5,4))
  
  # This just copies the function as written, so not really an independent check.
  #  Should pass by definition. But will help in case we change the correlation
  #  function later, or convert it to Rcpp.
  wendland2func <- function(x1,x2,theta) {
    diffmat =abs(outer(x1,x2,'-'))
    expLS = exp(3*theta[1])
    h = diffmat/expLS
    C = pmax(1 - h, 0)^5 * (8*h^2 + 5*h + 1)
    C
  }
  corr2 <- wendland2func(x1, x2, theta=th)
  expect_equal(corr1, corr2)
  
  # Now check that dC is actually grad of C
  corr_C_dC <- CGGP_internal_CorrMatWendland2(x1=x1, x2=x2, theta=th, return_dCdtheta=TRUE)
  eps <- 1e-6
  for (i in 1:1) {
    thd <- c(0)
    thd[i] <- eps
    numdC <- (CGGP_internal_CorrMatWendland2(x1=x1, x2=x2, theta=th+thd) -
                CGGP_internal_CorrMatWendland2(x1=x1, x2=x2, theta=th-thd)) / (2*eps)
    # Should be more accurate but was exactly the same
    # plot(numdC, corr_C_dC$dCdtheta[,(1+4*i-4):(4*i)])
    expect_equal(numdC, corr_C_dC$dCdtheta[,(1+4*i-4):(4*i)], info = paste("theta dimension with error is",i))
  }
})


test_that("Logs work for all", {
  corrs <- list(CGGP_internal_CorrMatCauchySQ, CGGP_internal_CorrMatCauchySQT,
             CGGP_internal_CorrMatCauchy, CGGP_internal_CorrMatGaussian,
             CGGP_internal_CorrMatMatern32, CGGP_internal_CorrMatMatern52,
             CGGP_internal_CorrMatPowerExp, CGGP_internal_CorrMatWendland0,
             CGGP_internal_CorrMatWendland1, CGGP_internal_CorrMatWendland2
  )
  # Some corr funcs work better on log scale (like standard Gaussian and Matern),
  # but some work better on normal scale (Wendland). Wendland corr funcs 
  # were showing giving errors
  use_log_scales <- c(rep(T, 7),
                      rep(F, 3))
  # use_log_scales <- rep(T, 10)
  n1 <- 5
  n2 <- 6
  for (icorr in rev(1:length(corrs))) {
    corr <- corrs[[icorr]]
    use_log_scale <- use_log_scales[icorr]
    numpara <- corr(return_numpara = T)
    x1 <- runif(n1)
    x2 <- runif(n2)
    theta <- runif(numpara)*2-1
    
    # Check that it actually equals log of normal corr
    c1 <- corr(x1 = x1, x2 = x2, theta = theta)
    c1_log <- corr(x1 = x1, x2 = x2, theta = theta, returnlogs = T)
    expect_is(c1, "matrix")
    expect_is(c1_log, "matrix")
    expect_equal(log(c1), c1_log)
    
    # Check grad has correct C
    d1_log <- corr(x1 = x1, x2 = x2, theta = theta, returnlogs = T, return_dCdtheta = T)
    expect_is(d1_log$dCdtheta, "matrix")
    expect_is(d1_log, "list")
    expect_equal(d1_log$C, c1_log)
    
    # Check grad matches, this is repeat of what is done in each section above
    corr_C_dC <- corr(x1 = x1, x2 = x2, theta = theta, return_dCdtheta=T)
    eps <- 1e-6
    for (i in 1:numpara) {
      thd <- rep(0, numpara)
      thd[i] <- eps
      numdC <- (corr(x1=x1, x2=x2, theta=theta+thd) -
                  corr(x1=x1, x2=x2, theta=theta-thd)) / (2*eps)
      # Should be more accurate but was exactly the same
      expect_equal(numdC, corr_C_dC$dCdtheta[,(1+n2*i-n2):(n2*i)], 
                   info = paste("theta dimension with error is",i, ", icor is", icorr))
    }
    
    # Check grad matches on log scale
    corr_C_dC_logs <- corr(x1 = x1, x2 = x2, theta = theta, return_dCdtheta=T,
                           returnlogs=use_log_scale)
    eps <- 1e-6
    for (i in 1:numpara) {
      thd <- rep(0, numpara)
      thd[i] <- eps
      numdC <- (corr(x1=x1, x2=x2, theta=theta+thd, returnlogs = use_log_scale) -
                  corr(x1=x1, x2=x2, theta=theta-thd, returnlogs = use_log_scale)) / (2*eps)
      # Should be more accurate but was exactly the same
      # For Wendland, need to convert NaN to 0's
      numdC <- ifelse(is.nan(numdC), 0, numdC)
      # plot(numdC, corr_C_dC_logs$dCdtheta[,(1+n2*i-n2):(n2*i)])
      expect_equal(numdC, corr_C_dC_logs$dCdtheta[,(1+n2*i-n2):(n2*i)],
                   info = paste("theta dimension with error is", i, ", icor is", icorr))
    }
    
    rm(numpara, c1, c1_log)
  }
})
