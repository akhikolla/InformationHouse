library(activegp)
library(hetGP)
context("grad_expr")

covtypes <- c("Gaussian", "Matern3_2", "Matern5_2")

### There are effectively three cases for Wij calculation:
### computing i == j element, derivative wrt i
### computing i == j element, derivative wrt k != i
### computing i != j element, derivative wrt i or j
### computing i != j element, derivative wrt k != i or j
### These cases need to be verified for differentiation vs args 1 and 2.

for (ct in 1:length(covtypes)) {
  test_that(paste(covtypes[ct], "Kernel Univariate Derivatives Work wrt first arg"), {

            a <- 0.7# First point
            b <- 0.3# Second point
            l <- 1# Lengthscale
            h <- 1e-8# Step size for finite difference approx
            thresh <- 1e-3#Allowable gap between finite differencing and analytic soltn.

            ##### w_ii
            expect_true(abs((activegp:::w_ii_cpp(a+h,b,l,ct) - activegp:::w_ii_cpp(a,b,l,ct)) / h - 
                            activegp:::grad_w_ii_cppa(a,b,l,ct)) < thresh)

            ##### w_ij
            # For a > b
            expect_true(abs((activegp:::w_ij_cpp(a+h,b,l,ct) - activegp:::w_ij_cpp(a,b,l,ct)) / h - 
                            activegp:::grad_w_ij_cppa(a,b,l,ct)) < thresh)
            # For a <= b
            temp <- a
            a <- b
            b <- temp
            expect_true(abs((activegp:::w_ij_cpp(a+h,b,l,ct) - activegp:::w_ij_cpp(a,b,l,ct)) / h - 
                            activegp:::grad_w_ij_cppa(a,b,l,ct)) < thresh)
            temp <- a
            a <- b
            b <- temp

            ##### Ikk
            # For a > b
            expect_true(abs((activegp:::Ikk_cpp(a+h,b,l,ct) - activegp:::Ikk_cpp(a,b,l,ct)) / h - 
                            activegp:::grad_Ikk_cppa(a,b,l,ct)) < thresh)
            # For a <= b
            temp <- a
            a <- b
            b <- temp
            expect_true(abs((activegp:::Ikk_cpp(a+h,b,l,ct) - activegp:::Ikk_cpp(a,b,l,ct)) / h - 
                            activegp:::grad_Ikk_cppa(a,b,l,ct)) < thresh)
            temp <- a
            a <- b
            b <- temp
  })
  test_that(paste(covtypes[ct], "Kernel Univariate Derivatives Work wrt second arg"), {

            a <- 0.7# First point
            b <- 0.3# Second point
            l <- 1# Lengthscale
            h <- 1e-5# Step size for finite difference approx
            thresh <- 1e-3#Allowable gap between finite differencing and analytic soltn.

            ##### w_ii
            expect_true(abs((activegp:::w_ii_cpp(a,b+h,l,ct) - activegp:::w_ii_cpp(a,b,l,ct)) / h - 
                            activegp:::grad_w_ii_cppb(a,b,l,ct)) < thresh)

            ##### w_ij
            # For a > b
            expect_true(abs((activegp:::w_ij_cpp(a,b+h,l,ct) - activegp:::w_ij_cpp(a,b,l,ct)) / h - 
                            activegp:::grad_w_ij_cppb(a,b,l,ct)) < thresh)
            # For a <= b
            temp <- a
            a <- b
            b <- temp
            expect_true(abs((activegp:::w_ij_cpp(a,b+h,l,ct) - activegp:::w_ij_cpp(a,b,l,ct)) / h - 
                            activegp:::grad_w_ij_cppb(a,b,l,ct)) < thresh)
            temp <- a
            a <- b
            b <- temp

            ##### Ikk
            # For a > b
            expect_true(abs((activegp:::Ikk_cpp(a,b+h,l,ct) - activegp:::Ikk_cpp(a,b,l,ct)) / h - 
                            activegp:::grad_Ikk_cppb(a,b,l,ct)) < thresh)
            # For a <= b
            temp <- a
            a <- b
            b <- temp
            expect_true(abs((activegp:::Ikk_cpp(a,b+h,l,ct) - activegp:::Ikk_cpp(a,b,l,ct)) / h - 
                            activegp:::grad_Ikk_cppb(a,b,l,ct)) < thresh)
            temp <- a
            a <- b
            b <- temp
  })

  test_that(paste(covtypes[ct], " Wij Gradients work"), {
            set.seed(123)
            nvar <- 2
            n <- 10
            thresh <- 1e-3#Allowable gap between finite differencing and analytic soltn.

            design <- matrix(runif(n*nvar, 0, 1), nrow = n)
            response <- rnorm(n)

            model <- mleHomGP(design, response, lower = rep(1e-4, nvar), upper = rep(1,nvar),
                              noiseControl = list(g_bounds = c(1e-6, 1e-6)), covtype = covtypes[ct]) 
            C <- C_GP(model)

            xnew <- runif(nvar,0,1)
            xnew <- matrix(xnew, nrow = 1)

            ## Gradient wrt first argument, i == j
            i <- 1
            j <- 1
            theta <- rep(1, nvar)#Lengthscale vector.
            num <- c()
            for (k in 1:nvar) {
              h <- 1e-6#Finite difference step size

              xh <- rep(0,nvar)
              xh[k] <- h

              num <- rbind(num, ((drop(activegp:::W_kappa_ij2(xnew + xh, C$model$X0, theta = theta, i - 1, j - 1, ct = C$ct)) - 
                                  drop(activegp:::W_kappa_ij2(xnew, C$model$X0, theta = theta, i - 1, j - 1, ct = C$ct))) / h))
            } 

            expect_true(norm(num - activegp:::grad_W_kappa_ij2(xnew, C$model$X0, theta = theta, i - 1, j - 1, ct = C$ct)) < thresh)

            ## Gradient wrt first argument, i != j
            i <- 1
            j <- 2
            theta <- rep(1, nvar)#Lengthscale vector.
            num <- c()
            for (k in 1:nvar) {
              h <- 1e-6#Finite difference step size

              xh <- rep(0,nvar)
              xh[k] <- h

              num <- rbind(num, ((drop(activegp:::W_kappa_ij2(xnew + xh, C$model$X0, theta = theta, i - 1, j - 1, ct = C$ct)) - 
                                  drop(activegp:::W_kappa_ij2(xnew, C$model$X0, theta = theta, i - 1, j - 1, ct = C$ct))) / h))
            } 

            expect_true(norm(num - activegp:::grad_W_kappa_ij2(xnew, C$model$X0, theta = theta, i - 1, j - 1, ct = C$ct)) < thresh)

            ## Gradient wrt second argument, i==j
            i <- 1
            j <- 1
            theta <- rep(1, nvar)#Lengthscale vector.
            num <- c()
            for (k in 1:nvar) {
              h <- 1e-6#Finite difference step size

              xh <- rep(0,nvar)
              xh[k] <- h

              num <- rbind(num, ((drop(activegp:::W_kappa_ij2(C$model$X0, xnew + xh, theta = theta, i - 1, j - 1, ct = C$ct)) - 
                                  drop(activegp:::W_kappa_ij2(C$model$X0, xnew, theta = theta, i - 1, j - 1, ct = C$ct))) / h))
            } 

            expect_true(norm(num - activegp:::grad_W_kappa_ij2(xnew, C$model$X0, theta = theta, i - 1, j - 1, ct = C$ct)) < thresh)

            ## Gradient wrt second argument, i != j
            i <- 1
            j <- 2
            theta <- rep(1, nvar)#Lengthscale vector.
            num <- c()
            for (k in 1:nvar) {
              h <- 1e-6#Finite difference step size

              xh <- rep(0,nvar)
              xh[k] <- h

              num <- rbind(num, ((drop(activegp:::W_kappa_ij2(C$model$X0, xnew + xh, theta = theta, i - 1, j - 1, ct = C$ct)) - 
                                  drop(activegp:::W_kappa_ij2(C$model$X0, xnew, theta = theta, i - 1, j - 1, ct = C$ct))) / h))
            } 

            expect_true(norm(num - activegp:::grad_W_kappa_ij2_w2(xnew, C$model$X0, theta = theta, i - 1, j - 1, ct = C$ct)) < thresh)
  })

  #test_that(paste(covtypes[ct], " beta/gamma gradients work "), {
  #          set.seed(123)
  #          thresh <- 1e-3#Allowable gap between finite differencing and analytic soltn.

  #          # Set up problem
  #          nvar <- 2
  #          n <- 10
  #          design <- matrix(runif(n*nvar, 0, 1), nrow = n)
  #          response <- apply(design, 1, sum)
  #          model <- mleHomGP(design, response, lower = rep(1e-4, nvar), upper = rep(1,nvar),
  #                            noiseControl = list(g_bounds = c(1e-6, 1e-6)), covtype = covtypes[ct]) 
  #          C <- C_GP(model)
  #          xnew <- runif(nvar,0,1)
  #          xnew <- matrix(xnew, nrow = 1)

  #          # Common Precomputations
  #          if(is.null(nrow(xnew))) xnew <- matrix(xnew, nrow = 1)
  #          nvar <- ncol(xnew)
  #          kn1 <- cov_gen(xnew, C$model$X0, theta = C$model$theta, type = C$model$covtype)
  #          theta <- sqrt(C$model$theta/2)
  #          new_lambda <- predict(object = C$model, x = xnew, nugs.only = TRUE)$nugs/C$model$nu2_hat
  #          vn <- drop(1 - kn1 %*% tcrossprod(C$model$Ki, kn1)) + new_lambda + C$model$eps
  #          Kikn <- tcrossprod(C$model$Ki, kn1)
  #          Kiyn <- C$model$Ki %*% C$model$Z0 # Ki yn

  #          ## Check beta gradient
  #          # i == j
  #          i <- 1
  #          j <- 1
  #          num <- rep(NA, nvar)
  #          for (k in 1:nvar) {
  #            h <- 1e-5
  #            xh <- rep(0, nvar)
  #            xh[k] <- h
  #            num[k] <- (activegp:::get_betagamma(C, xnew + xh, i, j, kn1, theta, Kikn, Kiyn, vn)$beta - 
  #                       activegp:::get_betagamma(C, xnew, i, j, kn1, theta, Kikn, Kiyn, vn)$beta) / h
  #          }

  #          expect_true(norm(num - activegp:::get_betagamma(C, xnew, i, j, kn1, theta, Kikn, Kiyn, vn, grad = TRUE)$dbeta) < thresh)

  #          # i != j
  #          i <- 1
  #          j <- 2
  #          num <- rep(NA, nvar)
  #          for (k in 1:nvar) {
  #            h <- 1e-5
  #            xh <- rep(0, nvar)
  #            xh[k] <- h
  #            a <- (activegp:::get_betagamma(C, xnew + xh, i, j, kn1, theta, Kikn, Kiyn, vn) - activegp:::get_betagamma(C, xnew, i, j, kn1, theta, Kikn, Kiyn, vn)) / h
  #            num[k] <- a[1]
  #          }

  #          expect_true(norm(num - activegp:::get_betagamma(C, xnew, i, j, kn1, theta, Kikn, Kiyn, vn, grad = TRUE)[[1]]) < thresh)

  #          ## Check gamma gradient
  #          thresh <- 1e-2
  #          # i == j
  #          i <- 1
  #          j <- 1
  #          num <- rep(NA, nvar)
  #          for (k in 1:nvar) {
  #            h <- 1e-5
  #            xh <- rep(0, nvar)
  #            xh[k] <- h
  #            a <- (activegp:::get_betagamma(C, xnew + xh, i, j, kn1, theta, Kikn, Kiyn, vn) - activegp:::get_betagamma(C, xnew, i, j, kn1, theta, Kikn, Kiyn, vn)) / h
  #            num[k] <- a[2]
  #          }

  #          expect_true(norm(num - activegp:::get_betagamma(C, xnew, i, j, kn1, theta, Kikn, Kiyn, vn, grad = TRUE)[[2]]) < thresh)

  #          # i != j
  #          i <- 1
  #          j <- 2
  #          num <- rep(NA, nvar)
  #          for (k in 1:nvar) {
  #            h <- 1e-5
  #            xh <- rep(0, nvar)
  #            xh[k] <- h
  #            a <- (activegp:::get_betagamma(C, xnew + xh, i, j, kn1, theta, Kikn, Kiyn, vn) - activegp:::get_betagamma(C, xnew, i, j, kn1, theta, Kikn, Kiyn, vn)) / h
  #            num[k] <- a[2]
  #          }

  #          expect_true(norm(num - activegp:::get_betagamma(C, xnew, i, j, kn1, theta, Kikn, Kiyn, vn, grad = TRUE)[[2]]) < thresh)
  #})

  test_that(paste(covtypes[ct], "C var gradients work"), {
            set.seed(123)
            nvar <- 3
            n <- 10
            thresh <- 1e-3#Allowable gap between finite differencing and analytic soltn.

            design <- matrix(runif(n*nvar, 0, 1), nrow = n)
            response <- apply(design, 1, sum)

            model <- mleHomGP(design, response, lower = rep(1e-4, nvar), upper = rep(1,nvar),
                              noiseControl = list(g_bounds = c(1e-6, 1e-6)), covtype = covtypes[ct]) 
            C <- C_GP(model)

            xnew <- runif(nvar,0,1)
            xnew <- matrix(xnew, nrow = 1)

            ## Check gradient
            num <- rep(NA, nvar)
            for (k in 1:nvar) {
              h <- 1e-5
              xh <- rep(0, nvar)
              xh[k] <- h
              num[k] <- (C_var(C, xnew + xh) - C_var(C, xnew)) / h
            }
            expect_true(sum((num - C_var(C, xnew, grad = TRUE))^2) < thresh)
  })

  test_that(paste(covtypes[ct], "C var 2 gradients work"), {
            set.seed(123)
            nvar <- 3
            n <- 10
            thresh <- 1e-3#Allowable gap between finite differencing and analytic soltn.

            design <- matrix(runif(n*nvar, 0, 1), nrow = n)
            response <- apply(design, 1, sum)

            model <- mleHomGP(design, response, lower = rep(1e-4, nvar), upper = rep(1,nvar),
                              noiseControl = list(g_bounds = c(1e-6, 1e-6)), covtype = covtypes[ct]) 
            C <- C_GP(model)

            xnew <- runif(nvar,0,1)
            xnew <- matrix(xnew, nrow = 1)

            ## Check gradient
            num <- rep(NA, nvar)
            for (k in 1:nvar) {
              h <- 1e-5
              xh <- rep(0, nvar)
              xh[k] <- h
              num[k] <- (C_var2(C, xnew + xh) - C_var2(C, xnew)) / h
            }
            expect_true(sum((num - C_var2(C, xnew, grad = TRUE))^2) < thresh)
  })

  test_that(paste(covtypes[ct], "C Tr gradients work"), {
            set.seed(123)
            nvar <- 3
            n <- 10
            thresh <- 1e-3#Allowable gap between finite differencing and analytic soltn.

            design <- matrix(runif(n*nvar, 0, 1), nrow = n)
            response <- apply(design, 1, sum)

            model <- mleHomGP(design, response, lower = rep(1e-4, nvar), upper = rep(1,nvar),
                              noiseControl = list(g_bounds = c(1e-6, 1e-6)), covtype = covtypes[ct]) 
            C <- C_GP(model)

            xnew <- runif(nvar,0,1)
            xnew <- matrix(xnew, nrow = 1)

            ## Check gradient
            num <- rep(NA, nvar)
            for (k in 1:nvar) {
              h <- 1e-5
              xh <- rep(0, nvar)
              xh[k] <- h
              num[k] <- (C_tr(C, xnew + xh) - C_tr(C, xnew)) / h
            }
            expect_true(sum((num - C_tr(C, xnew, grad = TRUE))^2) < thresh)
  })
}

test_that("kernel_expr_stability", {
  # small t + same point
  expect_true(!any(is.nan(activegp:::grad_W_kappa_ij2(matrix(c(0.5,0.76), 1), matrix(c(0.5,0.76), 1), theta = c(0.00001, 0.001), i1 = 0, i2 = 0, ct = 1))))
  expect_true(!any(is.nan(activegp:::grad_W_kappa_ij2(matrix(c(0.5,0.76), 1), matrix(c(0.5,0.76), 1), theta = c(0.00001, 0.001), i1 = 0, i2 = 0, ct = 2))))
  expect_true(!any(is.nan(activegp:::grad_W_kappa_ij2(matrix(c(0.5,0.76), 1), matrix(c(0.5,0.76), 1), theta = c(0.00001, 0.001), i1 = 0, i2 = 0, ct = 3))))
  expect_true(!any(is.nan(activegp:::grad_W_kappa_ij2(matrix(c(0.5,0.76), 1), matrix(c(0.5,0.76), 1), theta = c(0.00001, 0.001), i1 = 0, i2 = 1, ct = 1))))
  expect_true(!any(is.nan(activegp:::grad_W_kappa_ij2(matrix(c(0.5,0.76), 1), matrix(c(0.5,0.76), 1), theta = c(0.00001, 0.001), i1 = 0, i2 = 1, ct = 2))))
  expect_true(!any(is.nan(activegp:::grad_W_kappa_ij2(matrix(c(0.5,0.76), 1), matrix(c(0.5,0.76), 1), theta = c(0.00001, 0.001), i1 = 0, i2 = 1, ct = 3))))
  
})


