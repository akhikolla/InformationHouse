loglik.4 <- function(pa,x,y,pam,paf) {
  nam <- 2*x[1]+x[2]
  naf <- 2*y[1]+y[2]
  nbm <- 2*x[3]+x[2]
  nbf <- 2*y[3]+y[2]
  logm <- sum(log(c(pam,1-pam,2))*c(nam,nbm,x[2])) 
  logf <- sum(log(c(paf,1-paf,2))*c(naf,nbf,y[2])) 
  loglik0 <- logm+logf
  nparam0 <- 2
  res <- c(loglik0,nparam0)
  return(res)
}
