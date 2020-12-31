fitmodel12 <-
function(x,y,pinit=NULL,tol=0.000001) {
  nAAm <- x[1]
  nABm <- x[2]
  nBBm <- x[3]
  
  nAAf <- y[1]
  nABf <- y[2]
  nBBf <- y[3]
  
  nAf <- 2*nAAf + nABf
  nAm <- 2*nAAm + nABm
  
  nA <- nAf + nAm
  
  nm <- sum(x)
  nf <- sum(y)
  
  n <- nm + nf
  
  if(!is.null(pinit))  pinit <- pinit else pinit <- nA/(2*n)
  
  p <- pinit
  
  h <- 1000
  i <- 1
  
  while(abs(h) > tol) {
    valfun <- fobj(p,nAAf,nABf,nBBf,nA,nm,nAf,n)
    valder <- fder(p,nAAf,nABf,nBBf,nA,nm,nAf,n)
    h <- -valfun/valder  
    p <- p + h
    i <- i + 1
  }
  
  fhat <- (nA - 2*p*n)/(2*p*nm- nA + nAf)
  
  return(list(p=p,fhat=fhat))
}
