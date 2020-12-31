library(activegp)
library(hetGP)
context("kernel_expr")

test_that("kernel_expr", {
  
  ## Initial R codes
  erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
  ## Check Gaussian kernel expressions
  w_ii <- function(a, b, t){
    1/(8*t^3) * ((2* (-2+a+b) * exp((-a^2 -b^2 -2 + 2 * a + 2 * b)/(2*t^2)) * t + exp(-(a-b)^2/(4 *t^2)) *  sqrt(pi) * ((a-b)^2-2* t^2) * erf((-2+a+b)/(2* t))) - (2* (a+b)* t * exp(-((a^2+b^2)/(2*t^2))) + exp(-(a-b)^2/(4* t^2)) * sqrt(pi)* ((a-b)^2-2* t^2)* erf((a+b)/(2* t))))
  }
  w_ij <- function(a, b, t){
    -((2 * (exp(-(a^2 + b^2)/(2*t^2)) - exp((-a^2 -b^2 + 2 *(a + b -1))/(2*t^2)))* t + (a-b) * exp(-(a-b)^2/(4 * t^2)) * sqrt(pi) * (erf((-2+a+b)/(2 * t)) - erf((a+b)/(2* t)))))/(4 * t)
  }
  IkG <- function(x, theta){
    sqrt(pi/2) * theta * (erf(x/(sqrt(2) * theta)) - erf((x - 1)/(sqrt(2) * theta)))
  }
  Ikk <- function(a, b, theta){
    theta <- 2*theta^2
    sqrt(2*pi*theta)/4 * exp(-(a-b)^2/(2*theta))*(erf((2 - (a + b))/(sqrt(2*theta))) + erf((a + b)/(sqrt(2*theta))))
  }
  W_func_nd <- function(design, theta, i1, i2){
    if(is.null(dim(design))) design <- matrix(design, ncol = 1)
    n <- nrow(design)
    
    W <- matrix(NA, n, n)
    
    if(i1 == i2){
      for(i in 1:n){
        for(j in i:n){
          W[i, j] <- W[j,i] <- w_ii(design[i, i1], design[j, i1], theta[i1]) * prod(mapply(Ikk, a = design[i, -i1], b = design[j, -i1], t = theta[-i1]))
        }
      }
    }else{
      for(i in 1:n){
        for(j in 1:n){
          if(length(theta) > 2) tmp <- prod(mapply(Ikk, a = design[i, -c(i1,i2)], b = design[j, -c(i1, i2)], t = theta[-c(i1, i2)])) else tmp <- 1
          W[i, j] <- w_ij(design[i, i1], design[j, i1], theta[i1]) * w_ij(design[j, i2], design[i, i2], theta[i2]) * tmp
        }
      }
    }
    return(W)
  }
  # w_func_nd <- function(design, theta, i1){
  #   if(is.null(dim(design))) design <- matrix(design, nrow = 1)
  #   n <- nrow(design)
  #   wvec <- rep(0, n)
  #   for(i in 1:n){
  #     if(length(theta) == 1) tmp <- 1 else tmp <- prod(mapply(IkG, theta = theta[-i1], x = design[i,-i1]))
  #     wvec[i] <- (cov_gen(matrix(design[i, i1]), matrix(1), theta = 2*theta[i1]^2) - cov_gen(matrix(design[i, i1]), matrix(0), theta = 2*theta[i1]^2)) * tmp
  #   } 
  #   return(wvec)
  # }
  
  set.seed(4)
  nvar <- 4
  n <- 80
  design <- matrix(runif(n * nvar), n, nvar)
  theta <- runif(nvar)
  
  # for(i in 1:nvar){
  #   expect_equal(w_func_nd(design = design, theta = theta, i1 = i),
  #                activegp::w_kappa_i(w = rep(0, n), design = design, theta = theta, i1 = i - 1, start = 0),
  #                tol = 1e-8)
  # }
  
  for(i in 1:nvar){
    for(j in 1:nvar){
      expect_equal(W_func_nd(design = design, theta = theta, i1 = i, i2 = j),
                   W_kappa_ij(design = design, theta = theta, i1 = i - 1, i2 = j - 1, ct = 1),
                   tol = 1e-8)
    }
  }
  
  ## Test kernel matrices agree
  expect_equal(W_kappa_ij(design = design, theta = theta, i1 = 1, i2 = 2, ct = 1)[3,],
               drop(activegp:::W_kappa_ij2(design1 = design[3,,drop = F], design2 = design, theta = theta, i1 = 1, i2 = 2, ct = 1)))
  
})



#################################################################
## Marginal Kernel Expressions
covtypes <- c("Gaussian", "Matern3_2", "Matern5_2")
ks <- list(
           function(x1, x2, l) exp(-(x1-x2)^2/(2*l^2)), # Squared Exponential
           function(x1, x2, l) (1 + sqrt(3) * abs(x1 - x2) / l) * exp(-sqrt(3) * abs(x1-x2)/l), # Matern 3/2
           function(x1, x2, l) (1 + sqrt(5) * abs(x1 - x2) / l + 5 * (x1-x2)^2 / (3*l^2)) * exp(-sqrt(5) * abs(x1-x2)/l) # Matern 5/2
           )
for (ct in 1:length(covtypes)) {
  test_that(paste(covtypes[ct], " kernel Expressions are marginally accurate"), {
            set.seed(123)
            ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
            ## Test marginal quantities
            l <- 0.1
            h <- 1e-6
            k <- ks[[ct]]
            covtype <- covtypes[ct]

            kd1 <- function(x1, x2, l) (k(x1+h, x2, l) - k(x1, x2, l)) / h
            kd2 <- function(x1, x2, l) (k(x1, x2+h, l) - k(x1, x2, l)) / h

            # Some example quatities.
            a <- 0.7
            b <- 0.3
            xs <- runif(1e5,0,1)

            thresh <- 1e-2

            # Verify each individual quanitity
            # The first quantity, integral of the product of two partials
            expect_true(abs(mean(sapply(xs, function(x) kd1(x, a, l) * kd2(b,x, l))) - 
              activegp:::w_ii_cpp(a, b, l, ct)) < thresh)
            expect_true(abs(mean(sapply(xs, function(x) kd1(x, b, l) * kd2(a,x, l))) - 
              activegp:::w_ii_cpp(b, a, l, ct)) < thresh)

            # The second quantity, integral of the product of a partial and the original kernel
            expect_true(abs(mean(sapply(xs, function(x) kd1(x, a, l) * k(b,x, l))) - 
              activegp:::w_ij_cpp(a, b, l, ct)) < thresh)
            expect_true(abs(mean(sapply(xs, function(x) kd1(x, b, l) * k(a,x, l))) - 
              activegp:::w_ij_cpp(b, a, l, ct)) < thresh)

            # The third quantity, integral of the product of two kernels
            expect_true(abs(mean(sapply(xs, function(x) k(x, a, l) * k(b,x, l))) - 
              activegp:::Ikk_cpp(a, b, l, ct)) < thresh)
            expect_true(abs(mean(sapply(xs, function(x) k(x, b, l) * k(a,x, l))) - 
              activegp:::Ikk_cpp(b, a, l, ct)) < thresh)
  })
}

#################################################################
## Join Kernel Expressions
covtypes <- c("Gaussian", "Matern3_2", "Matern5_2")
ks <- list(
           function(x1, x2, l) exp(-(x1-x2)^2/(2*l^2)), # Squared Exponential
           function(x1, x2, l) (1 + sqrt(3) * abs(x1 - x2) / l) * exp(-sqrt(3) * abs(x1-x2)/l), # Matern 3/2
           function(x1, x2, l) (1 + sqrt(5) * abs(x1 - x2) / l + 5 * (x1-x2)^2 / (3*l^2)) * exp(-sqrt(5) * abs(x1-x2)/l) # Matern 5/2
           )
for (ct in 1:length(covtypes)) {
  test_that(paste(covtypes[ct], " kernel Expressions are jointly accurate (in Wij calculation)"), {
            set.seed(123)
            ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
            ## Test Joint quantities
            k_marg <- ks[[ct]]
            k <- function(x1, x2, theta) sigma * prod(sapply(1:P, function(p) k_marg(x1[p], x2[p], theta[p])))
            kd1 <- function(x1, x2, dp, theta) {hv <- rep(0, P);hv[dp]<-h;(k(x1+hv, x2, theta) - k(x1, x2, theta)) / h}
            kd2 <- function(x1, x2, dp, theta) {hv <- rep(0, P);hv[dp]<-h;(k(x1, x2+hv, theta) - k(x1, x2, theta)) / h}
            kdd1 <- function(x1, x2, dp1, dp2, theta) {hv <- rep(0, P);hv[dp1]<-h;(kd1(x1, x2+hv, dp2, theta) - kd1(x1, x2, dp2, theta)) / h}

            N <- 20
            P <- 2
            theta <- c(1,2)
            design <- matrix(runif(N*P), nrow = N)

            h <- 1e-5
            sigma <- 0.6

            iters <- 1e3

            i <- 1
            j <- 1

            ## Get a numerical estimate via finite differencing and Monte Carlo
            outers <- list()
            for (iter in 1:iters) {
              x <- runif(P)
              kappa_i <- rep(NA, N)
              for (n in 1:N) {
                kappa_i[n] <- kd1(x, design[n,], i, theta)
              }
              kappa_j <- rep(NA, N)
              for (n in 1:N) {
                kappa_j[n] <- kd1(x, design[n,], j, theta)
              }

              outers[[iter]] <- kappa_i %*% t(kappa_j)
            }
            num_est <- Reduce(function(i,j) i + j, outers) / iters

            # Compare to our analytical result
            an_est <- W_kappa_ij(design, theta, i-1, j-1, ct)


            thresh <- 0.1
            expect_true(norm(num_est - 
              sigma^(2) * an_est) < thresh)


  })
}
