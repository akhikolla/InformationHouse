loglik.5 <- function(x,y,tracing=1) {
  a <- c(x,y)
  out <- hwe.modelE.sol(a,tracing=tracing)
  k <- out$k
  pamhat <- k[1]
  pafhat <- k[2]
  fall   <- k[3]
  pvecm0  <- c(pamhat^2 + pamhat*(1-pamhat)*fall,
              2*pamhat*(1-pamhat)*(1-fall),
              (1-pamhat)^2 + pamhat*(1-pamhat)*fall)
  ind.m <- !(x==0)
  logvecm0 <- log(pvecm0[ind.m])
  
  pvecf0 <- c(pafhat^2 + pafhat*(1-pafhat)*fall,
              2*pafhat*(1-pafhat)*(1-fall),
              (1-pafhat)^2 + pafhat*(1-pafhat)*fall)
  ind.f <- !(y==0)
  logvecf0 <- log(pvecf0[ind.f])
  
  loglik0 <- sum(x[ind.m]*logvecm0) + sum(y[ind.f]*logvecf0)
  nparam0 <- 3
  res <- c(loglik0,nparam0)
  return(res)
}
