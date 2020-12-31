#' Sample Predictor Distributions
#' 
#' Uses fitted BBL model to explore predictor distributions
#' 
#' After \code{bbl} fit, the resulting model is used by this function to sample 
#' predictor distributions in each response group and find the most likely
#' preditor set using MCMC.
#' 
#' @param object Object of class \code{bbl}
#' @param nsteps Total number of MC steps
#' @param verbose Verbosity level of output
#' @param progress.bar Display progress bar
#' @examples
#' titanic <- as.data.frame(Titanic)
#' b <- bbl(Survived~., data=titanic[,1:4], weights=titanic$Freq)
#' pxy <- mcSample(b)
#' pxy
#' @export
#' 
mcSample <- function(object, nsteps=1000, verbose=1, progress.bar=TRUE){
  
  if(!is(object,'bbl')) stop('Object is not of class bbl')
  xlevels <- object$xlevels
  nvar <- length(xlevels)
  groups <- object$groups
  ny <- length(groups)
  
  x <- matrix('', nrow=nvar, ncol=ny)
  rownames(x) <- names(xlevels)
  colnames(x) <- groups
  for(i in seq_len(nvar))
    x[i,] <- sample(xlevels[[i]], 1)  # random initial config.
  xmax <- x
  emax <- rep(-Inf, ny)
  names(emax) <- groups
  
  if(verbose>1) progress.bar <- FALSE
  if(progress.bar)
    pb <- txtProgressBar(style=3)
  for(istep in seq_len(nsteps)){
    for(ix in seq_len(nvar)){
      xi <- xlevels[[ix]]
      for(iy in seq_len(ny)){
        h <- object$coefficients$h[[iy]][[ix]]
        J <- object$coefficients$J[[iy]][[ix]]
        pr <- rep(1, length(xi))
        names(pr) <- xi
        for(k in seq(2, length(xi))){
          E <- 0
          if(!xi[k] %in% names(h)) next()
          E <- E + h[xi[k]]
          for(j in seq_len(nvar)){
            if(j==ix) next()
            if(!x[j,iy] %in% colnames(J[[j]])) next()
            E <- E + J[[j]][xi[k],x[j,iy]]
          }
          pr[k] <- exp(E)
        }
        x[ix,iy] <- sample(xi, size=1, prob=pr)
      }
      if(verbose>1) cat('istep = ',istep,', ix = ',ix,'\n')
    }
    for(iy in seq_len(ny)){
      e <- energy(x=x[,iy], h=object$coefficients$h[[iy]],
                  J=object$coefficients$J[[iy]])
      if(e > emax[iy]){
        emax[iy] <- e
        xmax[,iy] <- x[,iy]
      }
    }
#   if(progress.bar & istep %% 100==0)
    if(progress.bar) 
      setTxtProgressBar(pb, value=istep/nsteps)
  }
  if(progress.bar) close(pb)
  
  return(list(xmax=xmax, emax=emax))
}

# energy function
energy <- function(x, h, J){
  
  E <- 0
  nvar <- length(x)
  for(i in seq_along(x)){
    if(!x[i] %in% names(h[[i]])) next()
    E <- h[[i]][x[i]]
    if(i==nvar) next()
    for(j in seq(i+1,nvar)){
      if(!x[j] %in% colnames(J[[i]][[j]])) next()
      E <- E + J[[i]][[j]][x[i],x[j]]
    }
  }
  return(E)
}