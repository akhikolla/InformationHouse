context("testoutputparamdim")

sel.methods_totest <- c("UCB", "TS", "MAP")
# sel.methods_totest <- c("MAP") # Much faster

# -----------------------------------------------
# Test different output parameter dimensions
# -----------------------------------------------

# Use 3 dim output. 3rd func is LinComb of first two, so PCA should have 2 dim
f1 <- function(x){x[1]+x[2]^2 + cos(x[3]^2*2*pi*4) - 3.3}
f2 <- function(x){x[1]^1.3+.4*sin(6*x[2])+10}
f3 <- function(x) {f1(x) + .3*f2(x)}
f4 <- function(x) {-1.1*f1(x) + .8*f2(x)}
f5 <- function(x) {.2*f1(x)}
f <- function(x) {
  if (is.matrix(x)) {return(t(apply(x, 1, f)))}
  c(f1(x), f2(x), f3(x), f4(x), f5(x))#, rep(f5(x),20))
}
d <- 3
outd <- 5
outd_pca <- 2

nsup <- 20
xsup <- matrix(runif(nsup*d), nsup, d)
ysup <- f(xsup)
eps.sup <- 1e-2 # Use difference accuracy for supp data preds


ntest <- 20
xtest <- matrix(runif(ntest*d), ntest, d)
ytest <- f(xtest)
eps.test <- 1e-2 # Use difference accuracy for testp data preds


test_that("2. MV output, NO PCA, 1opd", {
  
  
  # First check MV with PCA
  SG <- CGGPcreate(d=d, batchsize=30)
  expect_is(SG, "CGGP")
  y <- f(SG$design)
  expect_error(SG <- CGGPfit(SG, Y=y), NA) # No error
  expect_length(SG$thetaMAP, d*SG$numpara)
  expect_true(!is.matrix(SG$thetaMAP))
  expect_true(ncol(SG$Y) == outd)
  expect_true(ncol(SG$y) == outd)
  yMVpred <- CGGPpred(SG, SG$design)$mean
  expect_equal(yMVpred, y, 1e-4)
  expect_equal(dim(yMVpred), c(nrow(SG$design), outd))
  
  # Check that append works without error, don't save it
  for (sel.method in sel.methods_totest) {
    expect_error(CGGPappend(SG, 30, sel.method), NA)
  }
  
  # Add supplemental data
  expect_error(SG <- CGGPfit(SG, Y=y, Xs=xsup, Ys=ysup), NA) # No error
  ysuppred <- CGGPpred(SG, xsup)$me
  expect_equal(ysuppred, ysup, eps.sup)
  
  # Check that append works with grid+supp data
  for (sel.method in sel.methods_totest) {
    expect_error(CGGPappend(SG, 30, sel.method), NA)
  }
})


test_that("4. MV output, NO PCA, separate opd", {
  
  
  # First check MV with PCA
  SG <- CGGPcreate(d=3, batchsize=30)
  expect_is(SG, "CGGP")
  y <- f(SG$design)
  
  # Fit a model to only each individual output, will compare thetaMAP later
  seed <- sample(1:10000, 1)
  set.seed(seed)
  SG1 <- CGGPfit(SG, Y=y[,1])
  SG2 <- CGGPfit(SG, Y=y[,2])
  SG3 <- CGGPfit(SG, Y=y[,3])
  SG4 <- CGGPfit(SG, Y=y[,4])
  SG5 <- CGGPfit(SG, Y=y[,5])
  
  # Now fit all
  set.seed(seed)
  expect_error(SG <- CGGPfit(SG, Y=y, separateoutputparameterdimensions = T), NA) # No error
  expect_length(SG$thetaMAP, d*SG$numpara*outd)
  expect_true(is.matrix(SG$thetaMAP))
  expect_true(ncol(SG$thetaMAP) == outd)
  expect_true(ncol(SG$Y) == outd)
  expect_true(ncol(SG$y) == outd)
  
  # Check that thetaMAP for each dimension matches when model is only fit to that dimension.
  expect_equal(SG$thetaMAP[,1], SG1$thetaMAP)
  expect_equal(SG$thetaMAP[,2], SG2$thetaMAP)
  expect_equal(SG$thetaMAP[,3], SG3$thetaMAP)
  expect_equal(SG$thetaMAP[,4], SG4$thetaMAP)
  expect_equal(SG$thetaMAP[,5], SG5$thetaMAP)
  # Check that predictions on these match
  expect_equal(c(CGGPpred(SG, xtest)$me[,1]), c(CGGPpred(SG1, xtest)$me))
  expect_equal(c(CGGPpred(SG, xtest)$me[,2]), c(CGGPpred(SG2, xtest)$me))
  expect_equal(c(CGGPpred(SG, xtest)$me[,3]), c(CGGPpred(SG3, xtest)$me))
  expect_equal(c(CGGPpred(SG, xtest)$me[,4]), c(CGGPpred(SG4, xtest)$me))
  expect_equal(c(CGGPpred(SG, xtest)$me[,5]), c(CGGPpred(SG5, xtest)$me))
  
  # Now check predictions
  yMVpred <- CGGPpred(SG$design, CGGP=SG)$mean
  expect_equal(yMVpred, y, 1e-4)
  expect_equal(dim(yMVpred), c(nrow(SG$design), outd))
  
  # Check that append works without error, don't save it
  for (sel.method in sel.methods_totest) {
    expect_error(CGGPappend(SG, 30, sel.method), NA)
  }
  
  # Add supplemental data
  expect_error(SG <- CGGPfit(SG, Y=y, Xs=xsup, Ys=ysup), NA) # No error
  ysuppred <- CGGPpred(SG, xsup)$me
  expect_equal(ysuppred, ysup, eps.sup)
  
  # Check that append works with grid+supp data
  for (sel.method in sel.methods_totest) {
    expect_error(CGGPappend(SG, 30, sel.method), NA)
  }
})

