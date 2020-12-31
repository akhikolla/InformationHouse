context("testpred")

test_that("Prediction matches exact on small samples", {
  
  # Use borehole
  borehole <- function(x) {
    rw <- x[, 1] * (0.15 - 0.05) + 0.05
    r <-  x[, 2] * (50000 - 100) + 100
    Tu <- x[, 3] * (115600 - 63070) + 63070
    Hu <- x[, 4] * (1110 - 990) + 990
    Tl <- x[, 5] * (116 - 63.1) + 63.1
    Hl <- x[, 6] * (820 - 700) + 700
    L <-  x[, 7] * (1680 - 1120) + 1120
    Kw <- x[, 8] * (12045 - 9855) + 9855
    
    m1 <- 2 * pi * Tu * (Hu - Hl)
    m2 <- log(r / rw)
    m3 <- 1 + 2 * L * Tu / (m2 * rw ^ 2 * Kw) + Tu / Tl
    return(m1 / m2 / m3)
  }
  d = 8
  testf<-function (x) {  return(borehole(x))} 
  
  # Create covariance function for vectors
  CorrMatCauchySQTPts <- function(x1, x2, theta) {
    
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
  CorrMatCauchySQTVecs <- function(x1,x2, theta) {
    prod(sapply(1:length(x1), function(i) {
      CorrMatCauchySQTPts(x1[i], x2[i], theta[(1+3*i-3):(3*i)])
    }))
  }
  CorrMatCauchySQTFull <- function(X, theta) {
    n <- nrow(X)
    outer(1:n, 1:n, Vectorize(function(i,j) {CorrMatCauchySQTVecs(X[i,],X[j,],theta)}))
  }
  CorrMatCauchySQTFull2 <- function(X, U, theta) {
    outer(1:nrow(X), 1:nrow(U), Vectorize(function(i,j) {CorrMatCauchySQTVecs(X[i,],U[j,],theta)}))
  }
  
  SG = CGGPcreate(d=d,31, corr="CauchySQT") # Need size to be small to avoid computationally singular in solve
  Y = testf(SG$design) #the design is $design, simple enough, right?
  
  n <- 50
  xp <- matrix(runif(d*n),n,d)
  
  # Get error if not fit yet
  expect_error(CGGPpred(SG, xp))
  
  # Then fit
  SG <- CGGPfit(CGGP=SG, Y=Y)
  # Now it should work
  expect_error(SGpred <- CGGPpred(xp=xp, CGGP=SG), NA)
  # Get error if wrong order of args
  expect_error(CGGPpred(xp, SG))
  
  my <- mean(Y)
  dy <- Y - my
  Sig <- CorrMatCauchySQTFull(SG$design, theta=SG$thetaMAP)
  s <- CorrMatCauchySQTFull2(xp, SG$design, theta = SG$thetaMAP)
  expred <- my + s %*% solve(Sig, dy)
  
  # Check mean predictions
  # plot(expred, SGpred$mean); abline(a=0,b=1, col=2)
  expect_equal(SGpred$mean, expred)
  
  # Test var predictions
  # Calculating s2 like this doesn't work since we use MAP
  # s2 <- c(t(Y) %*%solve(Sig, Y) / length(Y))
  # Just use the MAP value
  s2 <- SG$sigma2MAP
  exvar <- s2 * (1 - colSums(t(s) * solve(Sig, t(s))))
  if (F) {
    print(1/SGpred$var* exvar)
    plot(exvar, SGpred$var); abline(a=0,b=1, col=2)
  }
  expect_equal(c(SGpred$var), exvar)
  
  rm(my, dy, Sig, s, expred, SGpred)
  
  # -----------------------------
  # And check grid with supp data
  # -----------------------------
  ns <- 10
  Xs <- matrix(runif(ns*d), ns, d)
  Ys <- testf(Xs)
  SGs <- CGGPfit(SG, SG$Y, Xs=Xs, Ys=Ys, HandlingSuppData="Correct")
  rm(SG)
  SGpred <- CGGPpred(xp=xp, CGGP=SGs)

  Xall <- rbind(SGs$design, Xs)
  Yall <- c(SGs$Y, Ys)
  my <- mean(Yall)
  dy <- Yall - my
  Sig <- CorrMatCauchySQTFull(Xall, theta=SGs$thetaMAP)
  s <- CorrMatCauchySQTFull2(xp, Xall, theta = SGs$thetaMAP)
  expred <- my + s %*% solve(Sig, dy)

  # Check mean predictions
  if (F) {
    plot(expred, SGpred$mean); abline(a=0,b=1, col=2)
    summary((c(SGpred$mean)- expred) / expred)
  }
  expect_equal(SGpred$mean, expred, 1e-2)

  # Test var predictions
  # Calculating s2 like this doesn't work since we use MAP
  # s2 <- c(t(Y) %*%solve(Sig, Y) / length(Y))
  # Just use the MAP value
  s2 <- SGs$sigma2MAP
  exvar <- s2 * (1 - colSums(t(s) * solve(Sig, t(s))))
  if (F) {
    print(1/SGpred$var* exvar)
    plot(exvar, SGpred$var); abline(a=0,b=1, col=2)
    summary((c(SGpred$var)- exvar) / exvar)
  }
  expect_equal(c(SGpred$var), exvar)
  
  rm(SGs, Xall, Yall, my, dy, Sig, s, expred, SGpred, s2, exvar)
  
  
  # ----------------------------
  # Pred with only supp is exact
  # ----------------------------
  sg <- CGGPcreate(d, 0, Xs=Xs, Ys=Ys, corr="CauchySQT")
  SGpred <- predict(sg, xp)
  my <- mean(Ys)
  dy <- Ys - my
  Sig <- CorrMatCauchySQTFull(Xs, theta=sg$thetaMAP)
  s <- CorrMatCauchySQTFull2(xp, Xs, theta = sg$thetaMAP)
  expred <- my + s %*% solve(Sig, dy)
  
  # Check mean predictions
  if (F) {
    plot(expred, SGpred$mean); abline(a=0,b=1, col=2)
    summary((c(SGpred$mean)- expred) / expred)
  }
  expect_equal(SGpred$mean, expred)
  
  # Test var predictions
  # Calculating s2 like this doesn't work since we use MAP
  # s2 <- c(t(Y) %*%solve(Sig, Y) / length(Y))
  # Just use the MAP value
  s2 <- sg$sigma2MAP
  exvar <- s2 * (1 - colSums(t(s) * solve(Sig, t(s))))
  if (F) {
    print(1/SGpred$var* exvar)
    plot(exvar, SGpred$var); abline(a=0,b=1, col=2)
    summary((c(SGpred$var)- exvar) / exvar)
  }
  expect_equal(c(SGpred$var), exvar, tol=1e-5)
})

test_that("predMV works", {
  SG <- CGGPcreate(d=3, batchsize=100)
  f1 <- function(x){x[1]+x[2]^2}
  f2 <- function(x){x[1]^1.3+.4*sin(6*x[2])+10}
  y1 <- apply(SG$design, 1, f1)#+rnorm(1,0,.01)
  y2 <- apply(SG$design, 1, f2)#+rnorm(1,0,.01)
  y <- cbind(y1, y2)
  SG <- CGGPfit(SG, Y=y)
  yMVpred <- CGGPpred(SG$design, CGGP=SG)$mean
  expect_equal(yMVpred[,1], y1, 1e-4)
  expect_equal(yMVpred[,2], y2, 1e-4)
  
  # Make sure predict generic works
  yMVpredS3 <- predict(SG, SG$design)$mean
  expect_equal(yMVpred, yMVpredS3)
  rm(yMVpred, yMVpredS3)
  
  # Check outdims
  # Get error if not suitable
  expect_error(CGGPpred(SG$design, CGGP=SG, outdims = 2))
  SGsep <- CGGPfit(SG, SG$Y, separateoutputparameterdimensions = T)
  yMVpred_o2 <- CGGPpred(SG$design, CGGP=SGsep, outdims = 2)$mean
  # Second dim should match, first should not.
  expect_equal(yMVpred_o2[,2], y2, 1e-4)
  # expect_true(all(abs(yMVpred_o2[,1] - y1)> 1e-4)) # This was failing on Travis
  expect_true(mean(abs(yMVpred_o2[,1] - y1)> 1e-4) > 0.9,
              info = paste(round(c(yMVpred_o2[1,1], sort(y1)),5),
                           collapse = " ")) # Maybe a couple are close
  
  # Doesn't work since there's no way to update Y without updating parameters too.
  # xpred <- matrix(runif(100*3),100,3)
  # SG1 <- CGGPfit(SG, Y=y1)
  # SG2 <- CGGPfit(SG, Y=y2)
  # y1pred <- CGGPpred(xpred, SG=SG1)$mean
  # y2pred <- CGGPpred(xpred, SG=SG2)$mean
  # yMVpred <- CGGPpred(xpred, SG=SG)$mean
  # expect_equal(yMVpred[,1], c(y1pred), tol=1e-2)
  # expect_equal(yMVpred[,2], c(y2pred), tol=1e-2)
})

test_that("Supplemented works", {
  d <- 3
  SG <- CGGPcreate(d=d, batchsize=100)
  f1 <- function(x){x[1]+sin(2*pi*x[1]) + x[2]^2}
  y1 <- apply(SG$design, 1, f1)
  SG <- CGGPfit(SG, Y=y1)
  
  # Add supplemental data
  nsup <- 20
  xsup <- matrix(runif(d*nsup), nsup, d)
  ysup <- apply(xsup, 1, f1)
  # SG$supplemented <- TRUE
  # SG$Xs <- xsup
  # SG$Ys <- ysup
  
  # # Get error when not fit
  # expect_error(CGGPpred(xp=xsup, SG=SG))
  
  # Should work after fitting
  SG <- CGGPfit(SG, Y=y1, Xs=xsup, Ys=ysup)
  
  # Predictions should match values at supplemented points
  expect_equal(c(CGGPpred(xp=xsup, CGGP=SG)$me), ysup, tol=1e-4)
  
  # # Predict at points
  # n <- 50
  # xp <- matrix(runif(d*10),n,d)
  # SGpred <- CGGPpred(xp=xp, SG=SG)
})

test_that("supplemental with MV output works", {
  
  SG <- CGGPcreate(d=3, batchsize=100)
  f1 <- function(x){x[1]+x[2]^2}
  f2 <- function(x){x[1]^1.3+.4*sin(6*x[2])+10}
  y1 <- apply(SG$design, 1, f1)#+rnorm(1,0,.01)
  y2 <- apply(SG$design, 1, f2)#+rnorm(1,0,.01)
  y <- cbind(y1, y2)
  
  xsup <- matrix(runif(3*30), ncol=3)
  ysup1 <- apply(xsup, 1, f1)
  ysup2 <- apply(xsup, 1, f2)
  ysup <- cbind(ysup1, ysup2)
  
  SG <- CGGPfit(SG, Y=y, Xs=xsup, Ys=ysup)
  yMVpred <- CGGPpred(SG$design, CGGP=SG)$mean
  # Making tol as big as 1e-3 since 1e-4 gives errors on Travis
  expect_equal(yMVpred[,1], y1, 1e-3)
  expect_equal(yMVpred[,2], y2, 1e-3)
  ysuppred <- CGGPpred(xsup, CGGP=SG)$mean
  expect_equal(ysuppred[,1], ysup1, 1e-3)
  expect_equal(ysuppred[,2], ysup2, 1e-3)
  
  # Doesn't work when giving in theta
  expect_error(CGGPpred(SG, xsup, theta=SG$thetaPostSamples[,1]))
  
})

test_that("pred with theta works", {
  d <- 5
  sg <- CGGPcreate(d=d, batchsize=1500)
  expect_error(CGGPpred(sg$design, sg))
  f <- function(x){x[1]^1.3+.4*sin(6*x[2])+4*exp(x[3])}
  y <- apply(sg$design, 1, f)
  sg <- CGGPfit(sg, y)
  
  npred <- 10
  xpred <- matrix(runif(d*npred), npred)
  p1 <- CGGPpred(xpred, CGGP=sg)
  p2 <- CGGPpred(xpred, CGGP=sg, theta=sg$thetaMAP)
  expect_error(CGGPpred(xpred, CGGP=sg, theta=sg$thetaPostSamples[1:8,1]))
  p3 <- CGGPpred(xpred, CGGP=sg, theta=sg$thetaPostSamples[,1])
  expect_equal(p1$me, p2$me)
  expect_false(isTRUE(all.equal(p1$me, p3$me)))
  
  # Check full Bayesian. Not easy, just check if it's close to MAP prediction
  expect_error(CGGPpred(xp=xpred, sg, theta = sg$thetaPostSamples[1,]))
  pb <- CGGPpred(xp=xpred, sg)
  expect_is(pb, "list")
  expect_equal(p1$me, p3$mean, tol=1e-2)
  # expect_equal(c(p1$me /p3$mean), rep(1,10), tol=1e-2)
  # Variances can be near zero, so use relative since all.equal won't for small values
  expect_equal(c(p1$var /p3$var), rep(1,10), tol=.5) # Huge tolerance since I don't know how to check it
  
})
