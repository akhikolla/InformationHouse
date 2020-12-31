context("Test imputation")

test_that("Imputation works", {
  d <- 3
  f <- function(x) {sin(x[1]) * x[2]}
  cg <- CGGPcreate(d, 100)
  
  # 1 NA
  Y <- apply(cg$design, 1, f)
  Y[90] <- NA
  expect_error(CGGPfit(cg, Y), NA)
  
  # 5 NA
  Y <- apply(cg$design, 1, f)
  Y[84:88] <- NA
  expect_error(CGGPfit(cg, Y), NA)
  
  # 9 NA, the 9 after base block
  Y <- apply(cg$design, 1, f)
  Y[2:10] <- NA
  expect_error(CGGPfit(cg, Y), NA)
  
  # Impute last 15. Check predictions work and aren't NA.
  Y <- apply(cg$design, 1, f)
  Y[2:10] <- NA
  expect_error(cg2 <- CGGPfit(cg, Y), NA)
  expect_error(cg2pred <- CGGPpred(cg2, cg2$design*.5+.2), NA)
  expect_true(all(!is.na(cg2pred$mean)))
  expect_true(all(!is.na(cg2pred$var)))
})


test_that("Imputation works on multidim output", {
  d <- 3
  f <- function(x) {t1 <- sin(x[1]) * x[2] + .3*exp(x[3]); c(t1, t1^1.1)}
  cg <- CGGPcreate(d, 100)
  
  # 1 NA
  Y <- t(apply(cg$design, 1, f))
  Y[90,1] <- NA
  expect_error(CGGPfit(cg, Y), NA)
  
  # 5 NA
  Y <- t(apply(cg$design, 1, f))
  Y[84:88, 1:2] <- NA
  expect_error(CGGPfit(cg, Y), NA)
  
  # Impute some. Check predictions work and aren't NA.
  Y <- t(apply(cg$design, 1, f))
  Y[82:90,1] <- NA
  expect_error(cg2 <- CGGPfit(cg, Y), NA)
  expect_error(cg2pred <- CGGPpred(cg2, cg2$design*.5+.2), NA)
  expect_true(all(!is.na(cg2pred$mean)))
  expect_true(all(!is.na(cg2pred$var)))
})