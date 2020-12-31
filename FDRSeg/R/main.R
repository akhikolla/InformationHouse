# SMUCE
smuce <- function(Y, q, alpha = 0.1, r = round(50/min(alpha, 1-alpha)), sd = stepR::sdrobnorm(Y))
{
  if (missing(q)) {
    q <- simulQuantile(1-alpha, length(Y), r, "smuce")
  }
  ret   <- .smuce_cpp(Y, q, sd)
  ret
}

# FDRSeg
fdrseg <- function(Y, q, alpha = 0.1, r = round(50/min(alpha, 1-alpha)), sd = stepR::sdrobnorm(Y))
{
  if (missing(q)) {
    q <- simulQuantile(1-alpha, length(Y), r, "fdrseg");
  }
  qm    <- cummax(q)
  ret   <- .fdrseg_cpp(Y, q, qm, sd)
  ret
}

# D-FDRSeg
dfdrseg <- function(Y, q, alpha = 0.1, r = round(50/min(alpha, 1-alpha)),
                       convKern, sd = stepR::sdrobnorm(Y, lag=length(convKern)+1)) {
  if (missing(convKern)) {
    stop("D-FDRSeg: 'convKern' has to be specified!")
  }
  if (missing(q)) {
    q <- simulQuantile(1-alpha, length(Y), r, "dfdrseg", convKern)
  }
  qm    <- cummax(q)
  ret   <- .dfdrseg_cpp(Y, q, qm, sd, length(convKern))
  ret
}