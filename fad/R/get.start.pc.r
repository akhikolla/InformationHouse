
get.start.pc <- function(X,iSD,mu,q,lower=0.005)
{
  n <- nrow(X)
  # f  <- function(x,args) {
  #   z <- args[[2]]*x;
  #   as.numeric({args[[1]] %*% z - sum(args[[3]]*z)}/args[[4]])
  # }
  # ftrans <- function(x,args) {
  #   as.numeric(args[[2]]*{-args[[3]]*sum(x) + crossprod(args[[1]],x)}/args[[4]])
  # }
  # svdsres <- svds(A = f,Atrans=ftrans,k=q,dim=dim(X),args=list(X,iSD,mu,sqrt(n)))
  svdsres  <- svds(X,k=q,nu=0,nv=q,opts=list(center=mu,scale=1/iSD)) #svds_XmD(X,mu,iSD,q)

  L <- .postmdiag(svdsres$v, svdsres$d) ; #diag(svdsres$d, q, q)
  start <- 1 - rowSums(L^2)
  start[start <=0 ] <- lower;
  return(start)
}
