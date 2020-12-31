context("testsupponly")

sel.methods_totest <- c("MAP", "UCB", "TS")
# sel.methods_totest <- c("MAP") # Only MAP to save time

test_that("1. Create, append, predict with only supp, scalar out", {
  d <- 3
  # f <- function(x){1*(cos(x[3]*2*pi*2)*x[1]^1.2 + (1-x[1]^.8)*sin(pi*x[2]^2) + log(x[3]+.2))}
  f <- function(x){1*(cos(x[1]*2*pi*2)*x[1]^1.2 + (1-x[1]^.8)*sin(pi*x[2]^2) + (1+x[2])*log(x[1]+.2))}
  
  set.seed(0)
  nsup <- 30
  xsup <- matrix(runif(nsup*d), nsup, d)
  ysup <- apply(xsup, 1, f)
  ntest <- 500
  xtest <- matrix(runif(ntest*d), ntest, d)
  ytest <- apply(xtest, 1, f)
  
  # Error if give in unname args to supp_args
  expect_error(CGGPcreate(d, 0, Xs=xsup, Ys=ysup, supp_args = list(12)))
  
  # Create with only supp
  expect_error(s1 <- CGGPcreate(d, 0, Xs=xsup, Ys=ysup, corr="CauchySQ"), NA)
  expect_true(is.null(s1[["design"]]))
  expect_equal(s1$uoCOUNT, 0)
  expect_equal(s1$poCOUNT, 1)
  expect_true(all(s1$uo==0))
  expect_equal(s1$po[1,], rep(1,d))
  expect_true(all(s1$po[-1,]==0))
  
  # Predict with only supp on supp points
  expect_error(p1sup <- CGGPpred(s1, xsup), NA)
  expect_equal(c(p1sup$mean), ysup)
  expect_true(all(p1sup$var < 1e-8))
  
  # Predict with only supp on test points
  expect_error(p1test <- CGGPpred(s1, xtest), NA)
  expect_true(cor(c(p1test$mean), ytest) > .8)
  expect_true(all(p1test$var > 1e-8))
  expect_true(all(p1test$var < var(ytest)))
  # Check errors
  expect_error(CGGPpred(s1, xtest, theta = c(1,s1$thetaMAP))) # wrong theta length
  expect_equal(CGGPpred(s1, xtest, theta = s1$thetaMAP),
              CGGPpred(s1, xtest)) # Giving in MAP is same as normal pred
  expect_error(CGGPpred(s1, xtest, outdims=2)) # Can't give in outdims for 1od
  # Same preds for tiny change in theta (have to recalculate supppw and Sti)
  expect_equal(CGGPpred(s1, xtest, theta=s1$thetaMAP*.99999999)$m,
               CGGPpred(s1, xtest)$m, tol=1e-6)
  expect_equal(CGGPpred(s1, xtest, theta=s1$thetaMAP*.99999999)$v,
               CGGPpred(s1, xtest)$v, tol=1e-6)
  
  # Append points, all three methods
  set.seed(0) # Fails on test, but never in console
  for (sel.method in sel.methods_totest) {
    expect_error(s1.app <- CGGPappend(s1, 100, sel.method), NA)
    expect_true(nrow(s1.app$design) > 90)
    expect_true(nrow(s1.app$design) < 100) # Can't get 100 since first block is size 1
    expect_equal(s1.app$design, s1.app$design_unevaluated)
    expect_true(all(s1.app$uo[1,] == 1)) # First block is initial
    expect_true(sum(s1.app$uo[2,]) == d+1) # 2nd block only has one 2
    # Make sure 3rd dim is least explored
    s1.app.colMeans <- colMeans(s1.app$uo[1:s1.app$uoCOUNT,])
    # print(c(sel.method,s1.app.colMeans))
    expect_true(s1.app.colMeans[1]+.1 > s1.app.colMeans[3], info = paste("s1s3",s1.app.colMeans, sel.method, collapse = " "))
    expect_true(s1.app.colMeans[2]+.1 > s1.app.colMeans[3], info = paste("s2s3",s1.app.colMeans, sel.method, collapse = " "))
  }
  set.seed(Sys.time())
  rm(s1, s1.app, s1.app.colMeans)
  
  # Again create with only supp, but use MCMC
  expect_error(s1 <- CGGPcreate(d, 0, Xs=xsup, Ys=ysup,
                                supp_args=list(numPostSamples=7)), NA)
  expect_true(is.null(s1[["design"]]))
  expect_equal(s1$uoCOUNT, 0)
  expect_equal(s1$poCOUNT, 1)
  expect_true(all(s1$uo==0))
  expect_equal(s1$po[1,], rep(1,d))
  expect_true(all(s1$po[-1,]==0))
  
  tp <- capture.output(print(s1))
  expect_is(tp, "character")
  expect_gt(length(tp), 6)
})




test_that("2. Create, append, predict with only supp, MVout, no PCA, yes sepOPD", {
  d <- 5
  # f <- function(x){1*(cos(x[3]*2*pi*2)*x[1]^1.2 + (1-x[1]^.8)*sin(pi*x[2]^2) + log(x[3]+.2))}
  f1 <- function(x){1*(cos(x[1]*2*pi*2)*x[1]^1.2 + (1-x[1]^.8)*sin(pi*x[2]^2) + (1+x[2])*log(x[1]+.2))}
  f2 <- function(x){x[1]*log(.8+x[4]) + x[4]^.3*sin(2*pi*x[2])}
  f3 <- function(x) {.3*f1(x) + 1.7*f2(x)}
  f <- function(x){
    if (is.matrix(x)) {cbind(apply(x, 1, f1), apply(x, 1, f2), apply(x, 1, f3))}
    else {c(f1(x), f2(x), f3(x))}
  }
  d_out <- 3
  d_outpca <- 2
  
  nsup <- 40
  xsup <- matrix(runif(nsup*d), nsup, d)
  ysup <- f(xsup)
  ntest <- 500
  xtest <- matrix(runif(ntest*d), ntest, d)
  ytest <- f(xtest)
  
  # Create with only supp
  expect_error(s1 <- CGGPcreate(d, 0, corr="Cauchy",
                                Xs=xsup, Ys=ysup, supp_args=list(separateoutputparameterdimensions=TRUE)), NA)
  expect_true(is.null(s1[["design"]]))
  expect_equal(s1$uoCOUNT, 0)
  expect_equal(s1$poCOUNT, 1)
  expect_true(all(s1$uo==0))
  expect_equal(s1$po[1,], rep(1,d))
  expect_true(all(s1$po[-1,]==0))
  # No PCA, no reduced dims
  expect_equal(ncol(s1$Ys), d_out)
  expect_equal(ncol(s1$ys), d_out)
  # theta correct dim
  expect_equal(ncol(s1$thetaMAP), d_out)
  
  # Predict with only supp on supp points
  expect_error(p1sup <- CGGPpred(s1, xsup), NA)
  expect_equal(p1sup$mean, ysup, tol=1e-6)
  # expect_true(all(p1sup$var < 1e-8))
  
  # Predict with only supp on test points
  expect_error(p1test <- CGGPpred(s1, xtest), NA)
  expect_true(all(diag(cor(p1test$mean, ytest)) > .8))
  # expect_true(all(p1test$var > 1e-8))
  # expect_true(all(p1test$var < var(ytest)))
  
  # Append points, all three methods
  for (sel.method in sel.methods_totest) {
    expect_error(s1.app <- CGGPappend(s1, 100, sel.method), NA)
    expect_true(nrow(s1.app$design) > 90)
    expect_true(nrow(s1.app$design) < 100) # Can't get 100 since first block is size 1
    expect_equal(s1.app$design, s1.app$design_unevaluated)
    expect_true(all(s1.app$uo[1,] == 1)) # First block is initial
    expect_true(sum(s1.app$uo[2,]) == d+1) # 2nd block only has one 2
    # Make sure 3rd dim is least explored
    s1.app.colMeans <- colMeans(s1.app$uo[1:s1.app$uoCOUNT,])
    expect_true(s1.app.colMeans[1]+.1 > s1.app.colMeans[3])
    expect_true(s1.app.colMeans[2]+.1 > s1.app.colMeans[3])
  }
  
  tp <- capture.output(print(s1))
  expect_is(tp, "character")
  expect_gt(length(tp), 6)
})


test_that("3. Create, append, predict with only supp, MVout, no PCA, no sepOPD", {
  d <- 5
  # f <- function(x){1*(cos(x[3]*2*pi*2)*x[1]^1.2 + (1-x[1]^.8)*sin(pi*x[2]^2) + log(x[3]+.2))}
  f1 <- function(x){1*(cos(x[1]*2*pi*2)*x[1]^1.2 + (1-x[1]^.8)*sin(pi*x[2]^2) + (1+x[2])*log(x[1]+.2))}
  f2 <- function(x){x[1]*log(.8+x[4]) + x[4]^.3*sin(2*pi*x[2])}
  f3 <- function(x) {.3*f1(x) + 1.7*f2(x)}
  f <- function(x){
    if (is.matrix(x)) {cbind(apply(x, 1, f1), apply(x, 1, f2), apply(x, 1, f3))}
    else {c(f1(x), f2(x), f3(x))}
  }
  d_out <- 3
  d_outpca <- 2
  nopd <- d_outpca
  
  nsup <- 40
  xsup <- matrix(runif(nsup*d), nsup, d)
  ysup <- f(xsup)
  ntest <- 500
  xtest <- matrix(runif(ntest*d), ntest, d)
  ytest <- f(xtest)
  
  # Create with only supp
  expect_error(s1 <- CGGPcreate(d, 0, corr="PowerExp",
                                Xs=xsup, Ys=ysup, supp_args=list(separateoutputparameterdimensions=FALSE)), NA)
  expect_true(is.null(s1[["design"]]))
  expect_equal(s1$uoCOUNT, 0)
  expect_equal(s1$poCOUNT, 1)
  expect_true(all(s1$uo==0))
  expect_equal(s1$po[1,], rep(1,d))
  expect_true(all(s1$po[-1,]==0))
  # No PCA, no reduced dims
  expect_equal(ncol(s1$Ys), d_out)
  expect_equal(ncol(s1$ys), d_out)
  # theta correct dim
  expect_is(s1$thetaMAP, "numeric")
  
  # Predict with only supp on supp points
  expect_error(p1sup <- CGGPpred(s1, xsup), NA)
  expect_equal(p1sup$mean, ysup, tol=1e-6)
  # expect_true(all(p1sup$var < 1e-8))
  
  # Predict with only supp on test points
  expect_error(p1test <- CGGPpred(s1, xtest), NA)
  expect_true(all(diag(cor(p1test$mean, ytest)) > .6)) # corr can go below .8
  # expect_true(all(p1test$var > 1e-8))
  # expect_true(all(p1test$var < var(ytest)))
  
  # Append points, all three methods
  for (sel.method in sel.methods_totest) {
    expect_error(s1.app <- CGGPappend(s1, 100, sel.method), NA)
    expect_true(nrow(s1.app$design) > 90)
    expect_true(nrow(s1.app$design) < 100) # Can't get 100 since first block is size 1
    expect_equal(s1.app$design, s1.app$design_unevaluated)
    expect_true(all(s1.app$uo[1,] == 1)) # First block is initial
    expect_true(sum(s1.app$uo[2,]) == d+1) # 2nd block only has one 2
    # Make sure 3rd dim is least explored
    s1.app.colMeans <- colMeans(s1.app$uo[1:s1.app$uoCOUNT,])
    expect_true(s1.app.colMeans[1]+.3 > s1.app.colMeans[3]) # Hard to get these 100%
    expect_true(s1.app.colMeans[2]+.3 > s1.app.colMeans[3])
  }
  
})

test_that("11. Single output: check that supp only matches grid predictions", {
  d <- 5
  f1 <- function(x){1*(cos(x[1]*2*pi*2)*x[1]^1.2 + (1-x[1]^.8)*sin(pi*x[2]^2) + (1+x[2])*log(x[1]+.2))}
  # f2 <- function(x){x[1]*log(.8+x[4]) + x[4]^.3*sin(2*pi*x[2])}
  # f3 <- function(x) {.3*f1(x) + 1.7*f2(x)}
  f <- function(x){
    if (is.matrix(x)) {apply(x, 1, f1)}
    else {f1(x)}
  }
  
  cg_grid <- CGGPcreate(d=d, 30)
  cg_grid <- CGGPfit(cg_grid, f(cg_grid$design))
  cg_supp <- CGGPcreate(d, 0, Xs=cg_grid$design, Ys=cg_grid$Y,
                        supp_args = list(set_thetaMAP_to=cg_grid$thetaMAP))
  expect_equal(cg_grid$thetaMAP, cg_supp$thetaMAP)
  
  xp <- matrix(runif(25*d), ncol=d)
  yp <- f(xp)
  if (F) {
    cbind(yp, predict(cg_grid, xp)$me, predict(cg_supp, xp)$me)
    cbind(yp, predict(cg_grid, xp)$var, predict(cg_supp, xp)$var)
    predict(cg_grid, xp)$me - predict(cg_supp, xp)$me
    plot(predict(cg_grid, xp)$me , predict(cg_supp, xp)$me); abline(a=0,b=1, col=2)
    plot(predict(cg_grid, xp)$va , predict(cg_supp, xp)$va); abline(a=0,b=1, col=2)
  }
  expect_equal(predict(cg_grid, xp)$me, predict(cg_supp, xp)$me)
  expect_equal(c(predict(cg_grid, xp)$va), predict(cg_supp, xp)$va, tol=.1)
  
  
  # Check likelihood matches. Won't actually match because of missing constants,
  #  but differences should match
  expect_equal(CGGP_internal_neglogpost(cg_supp$thetaMAP, cg_supp, Xs=cg_supp$Xs, ys=cg_supp$ys, y=NULL) - 
                 CGGP_internal_neglogpost(cg_supp$thetaMAP, cg_supp, Xs=cg_supp$Xs, ys=.5*cg_supp$ys^3, y=NULL),
               CGGP_internal_neglogpost(cg_grid$thetaMAP, cg_grid, cg_grid$y) - 
                 CGGP_internal_neglogpost(cg_grid$thetaMAP, cg_grid, .5*cg_grid$y^3)
  )
  expect_equal(CGGP_internal_gneglogpost(cg_supp$thetaMAP, cg_supp, Xs=cg_supp$Xs, ys=cg_supp$ys, y=NULL) - 
                 CGGP_internal_gneglogpost(cg_supp$thetaMAP, cg_supp, Xs=cg_supp$Xs, ys=.5*cg_supp$ys^3, y=NULL),
               c(CGGP_internal_gneglogpost(cg_grid$thetaMAP, cg_grid, cg_grid$y) - 
                   CGGP_internal_gneglogpost(cg_grid$thetaMAP, cg_grid, .5*cg_grid$y^3))
  )
})


test_that("12. MV, shared params, check that supp only matches grid predictions", {
  d <- 5
  f1 <- function(x){1*(cos(x[1]*2*pi*2)*x[1]^1.2 + (1-x[1]^.8)*sin(pi*x[2]^2) + (1+x[2])*log(x[1]+.2))}
  f2 <- function(x){x[1]*log(.8+x[4]) + x[4]^.3*sin(2*pi*x[2])}
  f3 <- function(x) {.3*f1(x) + 1.7*f2(x)}
  f <- function(x){
    if (is.matrix(x)) {cbind(apply(x, 1, f1), apply(x, 1, f2), apply(x, 1, f3))}
    else {c(f1(x), f2(x), f3(x))}
  }
  
  cg_grid <- CGGPcreate(d=d, 30)
  cg_grid <- CGGPfit(cg_grid, f(cg_grid$design), separateoutputparameterdimensions=F)
  cg_supp <- CGGPcreate(d, 0, Xs=cg_grid$design, Ys=cg_grid$Y,
                        supp_args = list(set_thetaMAP_to=cg_grid$thetaMAP,
                                         separateoutputparameterdimensions=F)
  )
  expect_equal(cg_grid$thetaMAP, cg_supp$thetaMAP)
  
  xp <- matrix(runif(25*d), ncol=d)
  yp <- f(xp)
  if (F) {
    cbind(yp, predict(cg_grid, xp)$me, predict(cg_supp, xp)$me)
    cbind(yp, predict(cg_grid, xp)$var, predict(cg_supp, xp)$var)
    predict(cg_grid, xp)$me - predict(cg_supp, xp)$me
    plot(predict(cg_grid, xp)$me, predict(cg_supp, xp)$me); abline(a=0,b=1,col=2)
    plot(predict(cg_grid, xp)$va, predict(cg_supp, xp)$va); abline(a=0,b=1,col=2)
  }
  expect_equal(predict(cg_grid, xp)$me, predict(cg_supp, xp)$me)
  expect_equal(predict(cg_grid, xp)$va, predict(cg_supp, xp)$va, tol=1e-1)
  
  # Check likelihood matches. Won't actually match because of missing constants,
  #  but differences should match
  expect_equal(CGGP_internal_neglogpost(cg_supp$thetaMAP, cg_supp, Xs=cg_supp$Xs, ys=cg_supp$ys, y=NULL) -
                 CGGP_internal_neglogpost(cg_supp$thetaMAP, cg_supp, Xs=cg_supp$Xs, ys=.5*cg_supp$ys^3, y=NULL),
               CGGP_internal_neglogpost(cg_grid$thetaMAP, cg_grid, cg_grid$y) -
                 CGGP_internal_neglogpost(cg_grid$thetaMAP, cg_grid, .5*cg_grid$y^3)
  )
  
  expect_equal(CGGP_internal_gneglogpost(cg_supp$thetaMAP, cg_supp, Xs=cg_supp$Xs, ys=cg_supp$ys, y=NULL) -
                 CGGP_internal_gneglogpost(cg_supp$thetaMAP, cg_supp, Xs=cg_supp$Xs, ys=.5*cg_supp$ys^3, y=NULL),
               CGGP_internal_gneglogpost(cg_grid$thetaMAP, cg_grid, cg_grid$y) -
                 CGGP_internal_gneglogpost(cg_grid$thetaMAP, cg_grid, .5*cg_grid$y^3)
  )
})


test_that("13. MV, separate params, check that supp only matches grid predictions", {
  d <- 5
  f1 <- function(x){1*(cos(x[1]*2*pi*2)*x[1]^1.2 + (1-x[1]^.8)*sin(pi*x[2]^2) + (1+x[2])*log(x[1]+.2))}
  f2 <- function(x){x[1]*log(.8+x[4]) + x[4]^.3*sin(2*pi*x[2])}
  f3 <- function(x) {.3*f1(x) + 1.7*f2(x)}
  f <- function(x){
    if (is.matrix(x)) {cbind(apply(x, 1, f1), apply(x, 1, f2), apply(x, 1, f3))}
    else {c(f1(x), f2(x), f3(x))}
  }
  
  cg_grid <- CGGPcreate(d=d, 30)
  cg_grid <- CGGPfit(cg_grid, f(cg_grid$design), separateoutputparameterdimensions = T)
  cg_supp <- CGGPcreate(d, 0, Xs=cg_grid$design, Ys=cg_grid$Y,
                        supp_args = list(set_thetaMAP_to=cg_grid$thetaMAP,
                                         separateoutputparameterdimensions=T)
  )
  expect_equal(cg_grid$thetaMAP, cg_supp$thetaMAP)
  
  xp <- matrix(runif(25*d), ncol=d)
  yp <- f(xp)
  if (F) {
    cbind(yp, predict(cg_grid, xp)$me, predict(cg_supp, xp)$me)
    cbind(yp, predict(cg_grid, xp)$var, predict(cg_supp, xp)$var)
    predict(cg_grid, xp)$me - predict(cg_supp, xp)$me
    plot(predict(cg_grid, xp)$me, predict(cg_supp, xp)$me); abline(a=0,b=1,col=2)
    plot(predict(cg_grid, xp)$va, predict(cg_supp, xp)$va); abline(a=0,b=1,col=2)
  }
  expect_equal(predict(cg_grid, xp)$me, predict(cg_supp, xp)$me)
  expect_equal(predict(cg_grid, xp)$va, predict(cg_supp, xp)$va, tol=1e-1)
  
  # Check likelihood matches. Won't actually match because of missing constants,
  #  but differences should match
  expect_equal(CGGP_internal_neglogpost(cg_supp$thetaMAP, cg_supp, Xs=cg_supp$Xs, ys=cg_supp$ys, y=NULL) -
                 CGGP_internal_neglogpost(cg_supp$thetaMAP, cg_supp, Xs=cg_supp$Xs, ys=.5*cg_supp$ys^3, y=NULL),
               CGGP_internal_neglogpost(cg_grid$thetaMAP, cg_grid, cg_grid$y) -
                 CGGP_internal_neglogpost(cg_grid$thetaMAP, cg_grid, .5*cg_grid$y^3),
               tol=1e-4
  )
})
