
fad.fit.cor <-
  function(R, q, start=NULL, maxit = 500, iSD=NULL, mu = NULL,lower = 0.005, control = NULL, ...)  {
    p <- ncol(R)
    if(is.null(start))
    {
      eigres <- eigs_sym(R,q)
      start <- 1 - rowSums( .postmdiag(eigres$vectors,sqrt(eigres$values))^2);
      start[start <= 0] = lower;
      start = cbind(start)
    }
    # print(start)
    fngr <- function(Psi) FAfngrC(Psi,R,q)

    res <- optim_fad(par = start,fngr = fngr,lower=lower,upper=1,method="L-BFGS-B",
                       control = c(list(fnscale=1,parscale = rep(0.01, length(start)),maxit = maxit,
                                        factr=1e2)))

    Lambda <- FAoutC(res$par, R, q,iSD,mu)
    dof <- 0.5 * ((p - q)^2 - p - q)
    un <- setNames(res$par, colnames(R))
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


FAfngrC <- function(Psi, R, q)
{
  p <- ncol(R)
  sc <- 1/sqrt(Psi)
  # # Args[[1]] = R; args[[2]] = sc
  # f  <- function(x,args) {
  #   z <- args[[2]]*x;
  #   v <- args[[1]] %*% z
  #   as.numeric(sc*v)
  # }
  if(!is.Matrix(R))
  {
    A <- sc*.postmdiag(R,sc);
  }  else
  {
    A <- sc*{R %*% Diagonal(x=sc)}
    A <- as(A,"symmetricMatrix")
  }
  eigsres <- eigs_sym(A = A,k=q)
  L <- eigsres$vectors

  e <- eigsres$values
  e <-  -sum(log(Psi)) - sum(1/Psi) -sum(log(e) - e) - q

  load <- .postmdiag(L , sqrt(pmax(eigsres$values-1,0)))
  load <- {1/sc}*load

  g <- rowSums(load^2) + Psi - 1;
  g <- g/Psi^2;
  return(list(-e,g));
}



FAoutC <- function(Psi, R, q,iSD,mu)
{
  p <- ncol(R)
  sc <- 1/sqrt(Psi)

  # Args[[1]] = R; args[[2]] = sc
  f  <- function(x,args) {
    z = args[[2]]*x;
    v = args[[1]] %*% z
    as.numeric(sc*v)
  }

  eigsres <- eigs_sym(A = f,k=q,n=p,args=list(R,sc))
  load <- .postmdiag(eigsres$vectors , sqrt(pmax(eigsres$values-1,0)))
  load <- {1/sc}*load

  return(load)
}


