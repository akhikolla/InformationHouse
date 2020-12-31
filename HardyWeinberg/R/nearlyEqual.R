nearlyEqual <- function(a,b,epsilon=1e-10) {
  if(length(a) != length(b)) stop("a and b have unequal length.")
  absA <- abs(a);
  absB <- abs(b);
  diff <- abs(a - b);
  y  <- rep(NA,length(a))
  i1 <- a == b
  i2 <- (a == 0 | b == 0 | (absA + absB < .Machine$double.xmin))
  y[i1] <- TRUE
  y[!i1 &  i2] <- (diff[!i1 & i2] < (epsilon*.Machine$double.xmin))
  y[!i1 & !i2] <- (diff[!i1 & !i2]/min((absA[!i1 & !i2] + absB[!i1 & !i2]), .Machine$double.xmax) < epsilon)
  return(y)
}
