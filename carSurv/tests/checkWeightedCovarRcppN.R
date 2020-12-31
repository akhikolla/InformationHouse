##########################################################
# Compare weighted unbiased variance C version with R code
library(carSurv)

# Standard R-Funktion
weightSelfCovar <- function(x, y, w) {
  meanX <- mean(x=x)
  wMeanY <- sum(w*y) / length(y)
  firstSum <- sum(w*(x-meanX)*(y-wMeanY))
  firstSum / length(x)
}

# Simple example: 
# Does Rcpp version have identical results with R-Version Function?
set.seed(1)
response <- rnorm(250)
set.seed(5)
response2 <- rnorm(250)
set.seed(2)
weightResp <- runif(249)
weightResp <- weightResp / sum(weightResp)
weightResp[250] <- 1-sum(weightResp)
check1 <- all.equal(weightedCovarRcppN(x=response, y=response2, w=weightResp),
                    weightSelfCovar(x=response, y=response2, w=weightResp))
stopifnot(check1)

# Is Rcpp version faster than the regular variant?
library(microbenchmark)
microBench1 <- microbenchmark(weightedCovarRcppN(x=response, y=response2, w=weightResp), 
                              weightSelfCovar(x=response, y=response2, w=weightResp), 
                              times=25)
check2 <- summary(microBench1)$expr [which.min(summary(microBench1)$mean)]==
  summary(microBench1)$expr[1]
stopifnot(check2)
