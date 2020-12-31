context("testplot")

test_that("Plots work", {
  
  SG <- CGGPcreate(d=3, batchsize=20, corr="m52")
  f <- function(x){x[1]+x[2]*log(x[1]+.1)+sin(2*pi*4*x[2]^2) + sqrt(x[1])*cos(2*pi*5*x[3])}
  y <- apply(SG$design, 1, f)
  # SG <- CGGPfit(SG, Y=y)
  expect_error(SG <- CGGPfit(SG, Y=y), NA)
  
  # Heat map
  pheat <- CGGPplotheat(SG)
  expect_is(pheat, "gg")
  expect_is(pheat, "ggplot")
  
  # Histogram
  phist <- CGGPplothist(SG)
  expect_is(phist, "gg")
  expect_is(phist, "ggplot")
  # # Don't actually want warning for log scale, but getting it, so expect it
  # expect_warning(print(phist))
  
  # Blockplot
  expect_error(bp <- CGGPplotblocks(SG), NA)
  expect_is(bp, "ggplot")
  # Blockplot from matrix
  mat <- matrix(c(1,1,1,2,2,1,2,2,1,3), ncol=2, byrow=TRUE)
  expect_error(bp <- CGGPplotblocks(mat), NA)
  expect_is(bp, "ggplot")
  # Error if not CGGP object or matrix
  expect_error(bp <- CGGPplotblocks(c(0,1,2)))
  # Single plot with 2D matrix
  
  # Validation plot
  Xval <- matrix(runif(3*100), ncol=3)
  Yval <- apply(Xval, 1, f)
  pval <- CGGPvalplot(CGGP=SG, Xval=Xval, Yval=Yval)
  expect_is(pval, "gg")
  expect_is(pval, "ggplot")
  
  # Val stats
  expect_error(vstats <- valstats(c(0,1,2), c(.01,.01,.01), c(0,1.1,1.9)), NA)
  expect_is(vstats, "data.frame")
  expect_equal(nrow(vstats), 1)
  rm(vstats)
  expect_error(vstats <- valstats(c(0,1,2), c(.01,.01,.01), c(0,1.1,1.9), MAE=TRUE), NA)
  rm(vstats)
  expect_error(vstats <- valstats(c(0,1,2), c(.01,.01,.01), c(0,1.1,1.9),
                                  metrics=function(a,b,c)mean(abs(a-c))),NA)
  rm(vstats)
  expect_error(vstats <- valstats(c(0,1,2), c(.01,.01,.01), c(0,1.1,1.9),
                                  metrics=list(function(a,b,c)mean(abs(a-c)))), NA)
  rm(vstats)
  expect_error(vstats <- valstats(c(0,1,2), c(.01,.01,.01), c(0,1.1,1.9),
                                  metrics=list(mae=function(a,b,c)mean(abs(a-c)))), NA)
  rm(vstats)
  expect_error(vstats <- valstats(predmean=c(0,1,2), Yval=c(0,1.1,1.9)), NA)
  expect_is(vstats, "data.frame")
  expect_equal(nrow(vstats), 1)
  rm(vstats)
  
  # CGGP Val stats
  vstats <- CGGPvalstats(SG, Xval, Yval)
  expect_is(vstats, "data.frame")
  expect_equal(nrow(vstats), 1)
  
  rm(Xval, Yval)
  # Multiple outputs
  SG2 <- CGGPcreate(d=3, batchsize=100)
  f1 <- function(x){x[1]+x[2]^2/(x[1]+1) +sin(x[1]*x[3])}
  f2 <- function(x){x[1]^1.3+.4*sin(2*pi*x[2])*x[1] + (x[1]+1)*exp(x[3])+10}
  y1 <- apply(SG2$design, 1, f1)#+rnorm(1,0,.01)
  y2 <- apply(SG2$design, 1, f2)#+rnorm(1,0,.01)
  y <- cbind(y1, y2)
  SG2 <- CGGPfit(SG2, Y=y)
  SG2sep <- CGGPfit(SG2, Y=y, separateoutputparameterdimensions = T)
  Xval2 <- matrix(runif(3*100), ncol=3)
  Yval2 <- cbind(apply(Xval2, 1, f1),
                apply(Xval2, 1, f2))
  # Error with Yval is bad
  expect_error(CGGPvalstats(SG2, Xval2, cbind(Yval2, Yval2)))
  sv <- CGGPvalstats(SG2, Xval2, Yval2)
  expect_is(sv, "data.frame")
  expect_equal(nrow(sv), 2)
  # expect_true(all(sv$RMSE<.1))
  expect_equal(nrow(CGGPvalstats(SG2, Xval2, Yval2, bydim=FALSE)), 1)
  expect_error(CGGPvalstats(SG2, Xval2, Yval2[,1]))
  # Val plot
  pval <- CGGPvalplot(CGGP=SG2, Xval=Xval2, Yval=Yval2)
  expect_is(pval, "gg")
  expect_is(pval, "ggplot")
  pval2 <- CGGPvalplot(CGGP=SG2, Xval=Xval2, Yval=Yval2, d=2)
  expect_is(pval2, "gg")
  expect_is(pval2, "ggplot")
  
  # Corr plot
  expect_error(p <- CGGPplotcorr(), NA)
  expect_is(p, "ggplot")
  expect_error(p <- CGGPplotcorr(CGGP_internal_CorrMatCauchySQ, theta=c(-.9,.8,.7,-.8)), NA)
  expect_is(p, "ggplot")
  expect_error(p <- CGGPplotcorr(SG), NA)
  expect_is(p, "ggplot")
  expect_error(p <- CGGPplotcorr(SG2), NA)
  expect_is(p, "ggplot")
  expect_error(p <- CGGPplotcorr(SG2sep), NA)
  expect_is(p, "ggplot")
  expect_error(p <- CGGPplotcorr(SG2sep, outdims = 2), NA)
  expect_is(p, "ggplot")
  rm(p)
  
  # slice plot
  # These should work fine, return ggplot
  expect_error(pp1 <- CGGPplotslice(SG), NA)
  expect_is(pp1, "ggplot")
  expect_error(pp1 <- CGGPplotslice(SG2), NA)
  expect_is(pp1, "ggplot")
  expect_error(pp1 <- CGGPplotslice(SG2, outdims = 2), NA)
  expect_is(pp1, "ggplot")
  expect_error(pp2 <- CGGPplotslice(SG, proj = c(0)), NA)
  expect_is(pp2, "ggplot")
  # Error if proj is not 1 or d dim
  expect_error(CGGPplotslice(SG, proj=c(1,1)))
  expect_error(CGGPplotslice(SG, proj=c(1,1,1,1)))
  # Error if not all evaluated
  expect_error(CGGPplotslice(CGGPappend(SG, 16)))
  rm(pp1, pp2)
  
  # Variogram plot
  expect_error(vario <- CGGPplotvariogram(SG), NA)
  expect_is(vario, "ggplot")
  rm(vario)
  expect_error(vario <- CGGPplotvariogram(SG2, facet = 1), NA)
  expect_is(vario, "ggplot")
  rm(vario)
  expect_error(vario <- CGGPplotvariogram(SG2, facet = 2), NA)
  expect_is(vario, "ggplot")
  rm(vario)
  expect_error(vario <- CGGPplotvariogram(SG2, facet = 3), NA)
  expect_is(vario, "ggplot")
  rm(vario)
  expect_error(vario <- CGGPplotvariogram(SG2sep, facet = 1), NA)
  expect_is(vario, "ggplot")
  rm(vario)
  expect_error(vario <- CGGPplotvariogram(SG2sep, facet = 1, outdims = 2), NA)
  expect_is(vario, "ggplot")
  rm(vario)
  expect_error(CGGPplotvariogram(SG2, facet = "not valid facet"))
  
  # Theta plot
  expect_error(tp <- CGGPplottheta(SG), NA)
  expect_is(tp, "ggplot")
  rm(tp)
  expect_error(tp <- CGGPplottheta(SG2), NA)
  expect_is(tp, "ggplot")
  rm(tp)
  expect_error(tp <- CGGPplottheta(SG2sep), NA)
  expect_is(tp, "ggplot")
  rm(tp)
  
  # Block selection plot
  expect_error(tp <- CGGPplotblockselection(SG), NA)
  expect_is(tp, "ggplot")
  rm(tp)
  expect_error(tp <- CGGPplotblockselection(SG, indims=1:2), NA)
  expect_is(tp, "ggplot")
  rm(tp)
  
  
  # neglogpost samples plot
  expect_error(tsamp <- CGGPplotsamplesneglogpost(SG), NA)
  expect_is(tsamp, "ggplot")
  rm(tsamp)
  expect_error(tsamp <- CGGPplotsamplesneglogpost(SG2), NA)
  expect_is(tsamp, "ggplot")
  rm(tsamp)
  expect_error(tsamp <- CGGPplotsamplesneglogpost(SG2sep), NA)
  expect_is(tsamp, "ggplot")
  rm(tsamp)
  # Force an Inf neglogpost
  SG$thetaPostSamples[,100] <- rep(10, SG$numpara)
  expect_warning(tsamp <- CGGPplotsamplesneglogpost(SG))
  expect_is(tsamp, "ggplot")
  rm(tsamp)
})
