
# evaluate step function
evalStepFun <- function(stepF)
{
  ret <- rep(0, stepF$n)
  left <- c(stepF$left, stepF$n)
  for (i in 1:length(stepF$left)) {
    ret[left[i]:left[i+1]] <- stepF$value[i]
  }
  ret
}

# create 'teeth' function
teethfun <- function(n, K, h=3) 
{
  u <- rep(0, n)
  s0 <- n/(K+1)
  i <- 1
  s <- s0
  while (i <= n) {
    u[i:round(s)] <- h*(1 - (round(s/s0)%%2))
    i <- round(s) + 1
    s <- s + s0
  }
  u
}


# false discovery proportion
computeFdp <- function(u, eJ) {
  eJ   <- eJ[(eJ <= length(u)) & (eJ > 1)]
  njmp <- length(eJ)
  if (njmp == 0) {
    fdr <- 0
  } else {
    d <- c(1, eJ, length(u)+1)
    nfd <- 0
    for (i in 1:njmp) {
      led <- min(ceiling((d[i]+d[i+1])/2), d[i+1]-1)
      red <- max(ceiling((d[i+1]+d[i+2])/2)-1, d[i+1])
      aux <- diff(u[led:red])
      if (all(aux == 0)) {
        nfd <- nfd + 1
      }
    }
    fdr <- nfd / (njmp+1)
  }
  fdr
}

# # borrow from "stepR" package
# # robust estimator of standard deviation for Gaussian data with jumps
# sdrobnorm <- function (x, p = c(0.25, 0.75), lag = 1)
# {
#   as.numeric( diff(quantile( diff(x, lag = lag), p)) / diff(qnorm(p)) / sqrt(2) )
# }

# V-measure, or validity-meansure
# Reference: 
#     [1] A. Rosenberg and J. Hirschberg, "V-Measure: A Conditional Entropy-Based External
#         Cluster Evaluation Measure.," EMNLP-CoNLL, no. June, pp. 410-420, 2007.
v_measure <- function(sig, est, beta=1)
{
  n <- length(sig)
  C <- c(0, which(diff(sig) != 0), n)
  K <- c(0, which(diff(est) != 0), n)
  
  # homogeneity
  HC <- .entropy_stepF(sig, n)
  if (HC == 0) {
    h <- 1
  } else {
    HCK <- 0
    for (i in 1:(length(K)-1)) {
      HCK <- HCK + .entropy_stepF(sig[(K[i]+1):K[i+1]], n)
    }
    h <- 1 - HCK / HC
  }
  
  # completeness
  HK <- .entropy_stepF(est, n)
  if (HK == 0) {
    c <- 1
  } else {
    HKC <- 0
    for (i in 1:(length(C)-1)) {
      HKC <- HKC + .entropy_stepF(est[(C[i]+1):C[i+1]], n)
    }
    c <- 1 - HKC / HK
  }
  
  # V_{\beta}
  (1 + beta)*h*c / (beta*h + c)
}

.entropy_stepF <- function(y, N) 
{
  n <- length(y)
  p <- diff(c(0, which(diff(y) != 0), n))
  sum(-p/N*log2(p/n))
} 
