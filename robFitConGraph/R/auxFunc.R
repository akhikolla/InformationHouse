myIPSC <- function(model,estimate,n=100,tol=10e-7) {
  diag(model) <- 0
  Laux <- myfitConGraphC(S=estimate,amat=model,n=n,tol=tol)
  df <- (sum(1-model) - ncol(model))/2
  SK <- estimate %*% solve(Laux[[1]])
  dev <- (sum(diag(SK)) - log(det(SK)) - ncol(model)) * n
  return(list(Shat=Laux[[1]], dev=dev, df=df, it=Laux[[2]]))
}

u <- function(s,p,d){return((d+p)/(d+s))}
