
fad.fit.X <-
  function(X, q, iSD=NULL, mu = NULL, start=NULL, lower = 0.005, control = NULL, maxit = 500, ...)
  {
    p <- ncol(X)
    n <- nrow(X)


    if(is.null(mu)) mu <- colMeans(X)
    if(is.null(iSD))
    {
      if(is.Matrix(X))
      {
        iSD <- 1/colSD(X,mu)/sqrt(n)
      } else      iSD <- apply(sweep(X,2,mu), 2,function(v) 1/sqrt(mean(v^2)))/sqrt(n)

    }

    if(is.null(start))
    {
      start <- get.start.pc(X,iSD,mu,q,lower=lower)
    }
    fngr <- function(Psi) FAfngr(Psi = Psi,X = X,q = q,iSD=iSD,mu=mu)

    res <- optim_fad(par = start,fngr = fngr,lower=lower,upper=1,method="L-BFGS-B",
                       control = c(list(fnscale=1,parscale = rep(0.01, length(start)),maxit = maxit,
                                        factr=1e2)))

    Lambda <- FAout(res$par, X, q,iSD,mu)
    dof <- 0.5 * ((p - q)^2 - p - q)
    un <- setNames(res$par, colnames(X))
    class(Lambda) <- "loadings"
    gerr = sqrt(sum({rowSums(Lambda^2)+res$par - 1}^2))
    ans <- list(converged = res$convergence == 0,
                loadings = Lambda, uniquenesses = un,
                gerr = gerr,
                criteria = c(objective = res$value, counts = res$counts),
                factors = q, dof = dof, method = "mle")
    class(ans) <- "fad"
    return(ans);
  }



# The following function assumes that Z'Z/n is a correlation matrix where
# Z = (X - mu)*diag(iSD)

FAfngr <- function(Psi, X, q,iSD,mu)
{
  p <- ncol(X)
  n <- nrow(X)
  sc <- sqrt(Psi)

  scaling <- sc/iSD #iSD/sc
  svdsres <- svds(X,k=q,nu=0,nv=q,opts=list(center=mu,scale=scaling)) # svds_XmD(X,mu,scaling,q)
  L <- svdsres$v

  e <- svdsres$d^2
  e <-  -sum(log(Psi)) - sum(1/Psi) -sum(log(e) - e) - q

  load <- .postmdiag(L , sqrt(pmax(svdsres$d^2-1,0))) #diag(sqrt(pmax(svdsres$d^2 - 1, 0)),q)
  load <- sc*load

  g <- rowSums(load^2) + Psi - 1;
  g <- g/Psi^2;
  return(list(-e,g));
}



FAout <- function(Psi, X, q,iSD,mu)
{
  p <- ncol(X)
  n <- nrow(X)
  sc <- sqrt(Psi)

  scaling <- sc/iSD #iSD/sc
  svdsres <- svds(X,k=q,nu=0,nv=q,opts=list(center=mu,scale=scaling)) # svds_XmD(X,mu,scaling,q)

  load <- .postmdiag(svdsres$v , sqrt(pmax(svdsres$d^2-1,0)))
  load <- sc*load
  return(load)
}

