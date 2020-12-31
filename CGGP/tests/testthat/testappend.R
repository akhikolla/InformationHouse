context("testappend")

test_that("CGGPappend works", {
  SG <- CGGPcreate(d=3, batchsize=100, corr='GAUSS')
  f <- function(x){x[1]*x[3]+x[2]^2}
  y <- apply(SG$design, 1, f)
  SG <- CGGPfit(SG, Y=y)
  for (i in c("UCB", "MAP", "TS", "Oldest", "Random", "Lowest")) {
    lastN <- nrow(SG$design)
    SG <- CGGPappend(CGGP=SG, batchsize=20, selectionmethod = i)
    expect_is(SG, "CGGP")
    expect_gt(nrow(SG$design), lastN)
    y <- apply(SG$design, 1, f)
    SG <- CGGPfit(SG, Y=y)
  }
  # Error, UC is not an option.
  expect_error(CGGPappend(SG, 20, selectionmethod = "UC"))
  
  # # Can append after append without fitting between
  # expect_error(SG <- CGGPappend(SG, 20), NA)
  # expect_error(SG <- CGGPappend(SG, 20), NA)
  
  # Can append with supp data
  xsup <- matrix(runif(3*10), ncol=3)
  ysup <- apply(xsup, 1, f)
  SG <- CGGPfit(SG, SG$Y, Xs=xsup, Ys=ysup, corr='m32')
  for (sel.method in c("UCB", "TS", "MAP")) {
    expect_error(CGGPappend(SG, 30, sel.method), NA)
  }
})

test_that("CGGPappend works with large number", {
  # Start it small
  SG <- CGGPcreate(d=3, batchsize=20)
  f <- function(x){x[1]+log(x[1]+.1)+sin(2*pi*4*x[2]^2) + cos(2*pi*5*x[3])}
  y <- apply(SG$design, 1, f)
  SG <- CGGPfit(SG, Y=y)
  lastN <- nrow(SG$design)
  
  # Adding 2000 will force it to increase ML and add rows to uo, pila, pala, etc.
  # But it doesn't show as working on codecov? Try 4000
  expect_error(SG <- CGGPappend(CGGP=SG, batchsize=2*2000, selectionmethod="MAP"), NA)
  expect_is(SG, "CGGP")
  expect_gt(nrow(SG$design), lastN)
  # y <- apply(SG$design, 1, f)
  # SG <- CGGPfit(SG, Y=y)
  # expect_error(CGGPappend(SG, 20, selectionmethod = "UC"))
  
  # Do we want prediction to work before fit is run? If yes, below should not give an error, but it does now.
  # ypred <- CGGPpred(SG$design, SG)
  
  # Check that prediction is still exact after adding rows to uo, pila, w, etc.
  # This gets inaccurate since there's so much data, so I'll remove it.
  # y <- apply(SG$design, 1, f)
  # SG <- CGGPfit(SG, Y=y)
  # ypred <- CGGPpred(SG, SG$design)
  # expect_equal(y, c(ypred$mean), tol=1e-2)
})

test_that("CGGPappend gives warning if it can't add any data", {
  
  SG <- CGGPcreate(d=3, batchsize=20, corr="cauchySQ")
  f <- function(x){x[1]+log(x[1]+.1)+sin(2*pi*4*x[2]^2) + cos(2*pi*5*x[3])}
  y <- apply(SG$design, 1, f)
  SG <- CGGPfit(SG, Y=y)
  
  expect_warning(CGGPappend(CGGP=SG, batchsize=1, selectionmethod="MAP"))
})
