loglik.3 <- function(pa,f,x,y,tracing=1) {
  a <- c(x,y)
#  fm <- ifelse(mac(x)==0,1,HWf(x))
 # ff <- ifelse(mac(y)==0,1,HWf(y))
#  ini <- c(pa,fm,ff)
  ini <- c(0.5,0,0)
  out <- hwe.modelC.sol(a,parinit=ini,tracing=tracing)
  k <- out$k
  pahat <- k[1]
  fm <- k[2]
  ff <- k[3]
  
  paam0 <- pahat^2 + pahat*(1-pahat)*fm
  pabm0 <- 2*pahat*(1-pahat)*(1-fm)
  pbbm0 <- (1-pahat)^2 + pahat*(1-pahat)*fm
  
  paaf0 <- pahat^2 + pahat*(1-pahat)*ff
  pabf0 <- 2*pahat*(1-pahat)*(1-ff)
  pbbf0 <- (1-pahat)^2 + pahat*(1-pahat)*ff
  
  vecm0 <- c(paam0,pabm0,pbbm0)
  vecf0 <- c(paaf0,pabf0,pbbf0)
  
  ind.x <- !(x==0)
  ind.y <- !(y==0)
  
  logvecm0 <- log(vecm0[ind.x])
  logvecf0 <- log(vecf0[ind.y])
  
  loglik0 <- sum(x[ind.x]*logvecm0) + sum(y[ind.y]*logvecf0)
  nparam0 <- 3
  res <- c(loglik0,nparam0)
  return(res)
}
