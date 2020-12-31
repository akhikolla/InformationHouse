# fastQuantile <- function (alpha, n, r=round(50/min(alpha, 1-alpha)), mType=c("norm-pen"), seed = 123, ...) {
#   set.seed(seed)
#   mType = match.arg(mType)
#   if (mType == "norm-pen") {
#     ms = numeric(r)
#     for (i in 1:r) {
#       cs    = cumsum(rnorm(n))
#       P     = data.frame("x1"=1:n, "x2"=cs)
#       Q     = data.frame("x1"=-(n-1):0, "x2"=c(-cs[(n-1):1],0))
#       Ri    = MinkowskiSum(P, Q, 1)
#       R     = P[Ri[,1],] + Q[Ri[,2],]
#       ms[i] = max(abs(R[,2])/sqrt(R[,1]) - sqrt(2 + 2*log(n/R[,1])))
#     }
#   } 
#   
#   quantile(ms, alpha,...)
# }

fastQuantile <- function(alpha, n, r=round(50/min(alpha, 1-alpha)), mType=c("norm-pen","pois"), seed = 123, ...) {
  ms = numeric(r)
  set.seed(seed)
  mType = match.arg(mType)
  if (mType == "norm-pen") {
    for (i in 1:r) {
      ms[i] = .mStat(rnorm(n), 0)
    }
  } else if (mType == "pois") {
    mu = 1
    for (i in 1:r) {
      ms[i] = .mStat(rpois(n,mu), 1)
    }
  } else {
    stop("Unsupported 'mType!")
  }
  
  quantile(ms, alpha,...)
}