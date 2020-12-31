# CGGPfit is used in many other tests, so not as many here

context("testfit")

test_that("CGGPfit works with Laplace approx", {
  d <- 3
  SG <- CGGPcreate(d=d, batchsize=30, corr = "m52")
  f <- function(x){x[1]+x[2]^2}
  y <- apply(SG$design, 1, f)
  # Check some fit input options
  expect_error(CGGPfit(SG, Y=as.matrix(y)), NA)
  expect_error(CGGPfit(SG, Ynew=as.matrix(y)), NA)
  expect_error(CGGPfit(SG, Y=c(y, 1.2)))
  # Set theta
  expect_error(CGGPfit(SG, Y=y, set_thetaMAP_to = SG$thetaMAP+.1), NA)
  # Give theta0
  expect_error(CGGPfit(SG, Y=y, theta0 = SG$thetaMAP+.1), NA)
  # Change correlation function
  expect_error(SG <- CGGPfit(SG, Y=y, corr = "cauchysqt"), NA)
  
  # Check that samples from Laplace are near max
  if (F) {
    stripchart(as.data.frame(t(SG$thetaPostSamples)))
    stripchart(as.data.frame(t(SG$thetaMAP)), add=T, col=2, pch=19)
  }
  expect_true(all(abs(rowMeans(SG$thetaPostSamples)- SG$thetaMAP) < apply(SG$thetaPostSamples, 1, sd)))
  
  # Check errors in neglogpost
  expect_error(CGGP_internal_neglogpost(rep(.5,length(SG$thetaMAP)), SG, y, HandlingSuppData="notanoption"))
  expect_error(CGGP_internal_neglogpost(rep(.5,length(SG$thetaMAP)), SG))
  expect_error(CGGP_internal_gneglogpost(rep(.5,length(SG$thetaMAP)), SG, y, HandlingSuppData="notanoption"))
  expect_error(CGGP_internal_gneglogpost(rep(.5,length(SG$thetaMAP)), SG))
  
  
  expect_error(neglogpost <- CGGP_internal_gneglogpost(rep(.5,length(SG$thetaMAP)), SG, y), NA)
  expect_is(neglogpost, "matrix")
  expect_length(neglogpost, length(SG$thetaMAP))
  
  
  expect_error(neglogpost2 <- CGGP_internal_gneglogpost(rep(.5,length(SG$thetaMAP)), SG, y, return_lik = T), NA)
  expect_is(neglogpost2, "list")
  expect_length(neglogpost2[[1]], 1)
  expect_length(neglogpost2[[2]], length(SG$thetaMAP))
  
  # Check neglogpost grad matches gneglogpost
  theta <- SG$thetaMAP / 2 # Don't want values near -1 or +1
  epsval <- 1e-4
  thetagrad <- CGGP_internal_gneglogpost(theta, SG, SG$y)
  for (i in 1:length(SG$thetaMAP)) {
    eps <- rep(0, length(SG$thetaMAP))
    eps[i] = eps[i] + epsval
    numgrad <- (CGGP_internal_neglogpost(theta + eps, SG, SG$y) - 
                  CGGP_internal_neglogpost(theta - eps, SG, SG$y)) / (2*epsval)
    expect_equal(thetagrad[i], numgrad, tol=1e-4)
  }
  
  if (F) {
    # expect_equal(
    #   c(CGGP_internal_gneglogpost(theta, SG, SG$y)),
    #   numDeriv::grad(CGGP_internal_neglogpost, theta, CGGP=SG, y=SG$y)
    # )
  }
  
  # Works with supplementary data
  nsup <- 30
  xsup <- matrix(runif(nsup*3), nsup, 3)
  ysup <- apply(xsup, 1, f)
  SG <- CGGPappend(SG, 20)
  ynew <- apply(SG$design_unevaluated, 1, f)
  expect_error(SG <- CGGPfit(SG, Ynew=ynew, Xs=xsup, Ys=ysup), NA)
  
  # Check neglogpost grad matches gneglogpost for all HandlingSuppData options
  theta <- SG$thetaMAP / 2 # Don't want values near -1 or +1
  epsval <- 1e-4
  for (handling in c("Correct", "Only", "Ignore")) {
    thetagrad <- CGGP_internal_gneglogpost(theta, SG, SG$y, Xs=xsup, ys=SG$ys, HandlingSuppData = handling)
    numgrad <- rep(0, length(SG$thetaMAP))
    for (i in 1:length(SG$thetaMAP)) {
      eps <- rep(0, length(SG$thetaMAP))
      eps[i] = eps[i] + epsval
      # numgrad <- (CGGP_internal_neglogpost(theta + eps, SG, SG$y, Xs=xsup, ys=SG$ys, HandlingSuppData = handling) - 
      #               CGGP_internal_neglogpost(theta - eps, SG, SG$y, Xs=xsup, ys=SG$ys, HandlingSuppData = handling)) / (2*epsval)
      numgrad[i] <- (-CGGP_internal_neglogpost(theta + 2*eps, SG, SG$y, Xs=xsup, ys=SG$ys, HandlingSuppData = handling) + 
                       8*CGGP_internal_neglogpost(theta + eps, SG, SG$y, Xs=xsup, ys=SG$ys, HandlingSuppData = handling) - 
                       8*CGGP_internal_neglogpost(theta - eps, SG, SG$y, Xs=xsup, ys=SG$ys, HandlingSuppData = handling) + 
                       CGGP_internal_neglogpost(theta - 2*eps, SG, SG$y, Xs=xsup, ys=SG$ys, HandlingSuppData = handling)) / (12*epsval)
      # 
      # print(numgrad)
    }
    expect_equal(c(thetagrad), numgrad, tol=1e-2, info = handling)
  }
  if (F) {
    # numDeriv::grad(function(th) {CGGP_internal_neglogpost(th, SG, SG$y, Xs=xsup, ys=SG$ys, HandlingSuppData = handling)}, theta)
  }
})

test_that("CGGPfit works with Ynew - scalar output", {
  
  
  SG <- CGGPcreate(d=3, batchsize=20)
  f <- function(x){x[1]+x[2]^2}
  y <- apply(SG$design, 1, f)
  # Can't give in Y and Ynew
  expect_error(SG <- CGGPfit(SG, Y=y, Ynew=y))
  # Works with just Ynew
  SG <- CGGPfit(SG, Ynew=y)
  expect_true(all(!is.na(SG$thetaPostSamples)))
  
  # After append, only give in Ynew
  SG <- CGGPappend(SG, 50)
  ynew <- apply(SG$design_unevaluated, 1, f)
  # Again error if give in Y and Ynew
  expect_error(CGGPfit(SG, Y=c(y, ynew), Ynew=ynew))
  # Get error when Ynew is wrong size
  expect_error(CGGPfit(SG, Ynew=ynew[1:(length(ynew)-1)]))
  # Works fine
  expect_is(SG <- CGGPfit(SG, Ynew = ynew), "CGGP")
})

test_that("CGGPfit works with Ynew - vector output", {
  
  
  SG <- CGGPcreate(d=3, batchsize=20)
  f1 <- function(x){x[1]+x[2]^2}
  f2 <- function(x){x[2]^1.3+sin(2*pi*x[3])}
  y1 <- apply(SG$design, 1, f1)
  y2 <- apply(SG$design, 1, f2)
  y <- cbind(y1, y2)
  # Can't give in Y and Ynew
  expect_error(SG <- CGGPfit(SG, Y=y, Ynew=y))
  # Other errors when giving ynew and wrong number of rows
  expect_error(CGGPfit(SG, Ynew=y[1:5,]))
  expect_error(CGGPfit(SG, Ynew=y1[1:5]))
  # Errors for bad theta0
  expect_error(CGGPfit(SG, Y=y, theta0 = matrix(0, ncol=3, nrow=SG$numpara), separateoutputparameterdimensions = T))
  expect_error(CGGPfit(SG, Y=y, theta0 = matrix(0, ncol=1, nrow=SG$numpara), separateoutputparameterdimensions = T))
  
  # Works with just Ynew
  SG <- CGGPfit(SG, Ynew=y)
  
  # After append, only give in Ynew
  SG <- CGGPappend(SG, 50)
  ynew1 <- apply(SG$design_unevaluated, 1, f1)
  ynew2 <- apply(SG$design_unevaluated, 1, f2)
  ynew <- cbind(ynew1, ynew2)
  # Again error if give in Y and Ynew
  expect_error(CGGPfit(SG, Y=rbind(y, ynew), Ynew=ynew))
  # Get error when Ynew is wrong size
  expect_error(CGGPfit(SG, Ynew=ynew[1:(nrow(ynew)-1),]))
  # Get error if you give in vector
  expect_error(CGGPfit(SG, Ynew=ynew1))
  # Works fine
  expect_is(SG <- CGGPfit(SG, Ynew = ynew), "CGGP")
  
  
})

# test_that("postvarmatcalc", {
#   # These don't work unless x1 and x2 have same length.
#   # Maybe only makes sense when x1 == x2?
#   x1 <- runif(5)
#   x2 <- x1
#   o1 <- CGGP_internal_postvarmatcalc(x1, x2,
#                                xo=c(.11), theta=c(.1,.2,.3),
#                                CorrMat=CGGP_internal_CorrMatCauchySQT,
#                                returndPVMC=F, returndiagonly=F)
#   o2 <- CGGP_internal_postvarmatcalc(x1, x2,
#                                      xo=c(.11), theta=c(.1,.2,.3),
#                                      CorrMat=CGGP_internal_CorrMatCauchySQT,
#                                      returndPVMC=F, returndiagonly=T)
#   o3 <- CGGP_internal_postvarmatcalc(x1, x2,
#                                      xo=c(.11), theta=c(.1,.2,.3),
#                                      CorrMat=CGGP_internal_CorrMatCauchySQT,
#                                      returndPVMC=T, returndiagonly=F)
#   o4 <- CGGP_internal_postvarmatcalc(x1, x2,
#                                      xo=c(.11), theta=c(.1,.2,.3),
#                                      CorrMat=CGGP_internal_CorrMatCauchySQT,
#                                      returndPVMC=T, returndiagonly=T)
#   expect_equal(o1, o3$Sigma_mat)
#   expect_equal(diag(o1), o2)
#   expect_equal(o4$Sigma_mat, diag(o3$Sigma_mat))
#   expect_equal(o4$dSigma_mat[,1], diag(o3$dSigma_mat[,1:5]))
#   expect_equal(o4$dSigma_mat[,2], diag(o3$dSigma_mat[,6:10]))
#   expect_equal(o4$dSigma_mat[,3], diag(o3$dSigma_mat[,11:15]))
# })

# test_that("", {
#   
# })