loglik.6 <- function(x,y,pam,paf,fm,ff) {
  ind.x <- !(x==0)
  ind.y <- !(y==0)
  pvec.m <- c(pam^2 + pam*(1-pam)*fm,
              2*pam*(1-pam)*(1-fm),
              (1-pam)^2 + pam*(1-pam)*fm)
  logvecm <- log(pvec.m[ind.x])
  
  pvec.f <- c(paf^2 + paf*(1-paf)*ff,
              2*paf*(1-paf)*(1-ff),
              (1-paf)^2 + paf*(1-paf)*ff)
  logvecf <- log(pvec.f[ind.y])
  loglik0 <- sum(x[ind.x]*logvecm) + 
             sum(y[ind.y]*logvecf)
  nparam0 <- 4
  res <- c(loglik0,nparam0)
  return(res)
}
