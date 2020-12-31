eh <- function(si, h, J){

  N <- length(si)
  e <- 0

  for(i in seq_len(N)){
    if(si[i]==0) next
    e <- e + h[[i]][si[i]]
    if(i < N){
      for(j in seq(i+1,N)){
        if(si[j]==0) next
        e <- e + J[[i]][[j]][si[i],si[j]]
      }
    }
  }
  return(exp(e))
}

enum <- function(si, L, i, h, J, e=NULL){

  for(s in seq(0,L[i]-1)){
    si[i] <- s
    if(i>1) e <- enum(si, L, i-1, h, J, e)
    else e <- rbind(e, cbind(t(si),eh(si, h, J)))
  }
  N <- length(si)
  return(e)
}

#' Generate Random Samples from Boltzmann Distribution
#' 
#' Random samples are drawn from Boltzmann distribution
#' 
#' All possible factor states are enumerated exhaustively using
#' input argument \code{predictors}. If the number of predictors \eqn{m}
#' or the number of factor levels \eqn{L_i} for each predictor \eqn{i}
#' are even moderately large (\eqn{m\ge 10} or \eqn{L_i\ge 5}), 
#' this function will likely hang because the number of all possible 
#' states grows exponentially.
#' 
#' @param nsample Sample size
#' @param predictors List of predictor factor levels.
#' @param h Bias parameter; see \code{\link{bbl}}.
#' @param J Interaction parameters; see \code{\link{bbl}}.
#' @param code_out Ouput in integer codes; \eqn{a_i = 0, \cdots, L_i-1}.
#'        If \code{FALSE}, output in factors in \code{predictors}.
#' @return Data frame of samples in rows and predictors in columns.
#' @examples
#' set.seed(512)
#' m <- 5
#' n <- 1000
#' predictors <- list()
#' for(i in 1:m) predictors[[i]] <- c('a','c','g','t')
#' par <- randompar(predictors)
#' xi <- sample_xi(nsample=n, predictors=predictors, h=par$h, J=par$J)
#' head(xi)
#' @export
sample_xi <- function(nsample=1, predictors=NULL, h, J, code_out=FALSE){

  L <- NULL
  for(p in predictors) L <- c(L, length(p))
  if(length(predictors)!=length(h) | length(predictors)!=length(J))
  stop("Predictors and h,J sizes don't match")
  nvar <- length(predictors)
  
  e <- enum(si=rep(0,nvar), L=L, i=nvar, h, J)
  sid <- sample(NROW(e), size=nsample, replace=TRUE, prob=e[,nvar+1])
  si <- e[sid,seq_len(nvar)]
  
  if(!code_out){
    for(i in seq_len(nvar)){
      x <- data.frame(x=factor(predictors[[i]][si[,i]+1], 
                                          levels=predictors[[i]]))
      if(i==1) fsi <- x
      else fsi <- cbind(fsi,x)
    }
  } else fsi <- as.data.frame(si)
  rownames(fsi) <- seq_len(nsample)
  colnames(fsi) <- names(predictors)
  return(fsi)
}

#' Generate Random Parameters
#' 
#' Random values of bias and interaction parameters are generated
#' using either uniform or normal distributions.
#' 
#' Input argument \code{predictors} is used to set up proper list 
#' structures of parameters.
#' 
#' @param predictors List of predictor factor levels. See \code{\link{bbl}}.
#' @param distr \code{c('unif','norm')} for uniform or normal distributions.
#' @param h0 Mean of bias parameters
#' @param dh \code{sd} of bias if \code{distr = 'unif'}. If \code{distr = 'norm'},
#'        \eqn{h = [h_0-dh, h_0+dh]}.
#' @param J0 Mean of interaction parameters.
#' @param dJ \code{sd} of interactions if \code{distr = 'unif'}. 
#'        If \code{distr = 'norm'}, \eqn{J = [J_0-dJ, J_0+dJ]}.
#' @return List of parameters, \code{h} and \code{J}.
#' @examples
#' set.seed(311)
#' predictors <- list()
#' for(i in 1:5) predictors[[i]] <- c('a','c')
#' par <- randompar(predictors=predictors)
#' par
#' @export
randompar <- function(predictors, distr='unif', h0=0, dh=1, J0=0, dJ=1){
  
  m <- length(predictors)
  L <- NULL
  for(p in predictors) L <- c(L, length(p))
  
  h <- J <- vector('list',m)
  for(i in seq_len(m)) 
    J[[i]] <- vector('list',m)
  
  for(i in seq_len(m)){
    if(distr=='unif')
      h[[i]] <- stats::runif(n=L[i]-1,min=h0-dh,max=h0+dh)
    else
      h[[i]] <- stats::rnorm(n=L[i]-1,mean=h0, sd=dh)
    names(h[[i]]) <- predictors[[i]][-1]
    for(j in seq(i,m)){
      if(i==j) x <- 0
      else{ 
        if(distr=='unif') 
          x <- stats::runif(n=(L[i]-1)*(L[j]-1), min=J0-dJ, max=J0+dJ)
        else x <- stats::rnorm(n=(L[i]-1)*(L[j]-1), mean=J0, sd=dJ)
      }
      x <- matrix(x, nrow=L[i]-1, ncol=L[j]-1)
      rownames(x) <- predictors[[i]][-1]
      colnames(x) <- predictors[[j]][-1]
      J[[i]][[j]] <- x
      if(i!=j) J[[j]][[i]] <- t(x)
    }
  }
  
  return(list(h=h, J=J))
}

#' Generate Random Boltzmann Bayes Model Data
#' 
#' Predictor-response paired data are generated
#' 
#' The argument \code{response} is used to set up all possible levels
#' of response groups and likewise for \code{predictors}. The parameter
#' argument \code{\link{par}} must have the appropriate structure 
#' consistent with \code{response} and \code{predictors}. This function
#' is a wrapper calling \code{\link{sample_xi}} multiple times.
#' 
#' @param predictors List of vectors of predictor levels
#' @param response Vector of response variables
#' @param prob Vector of probabilities for sampling each response group
#' @param par List of \code{\link{bbl}} parameters for each response group;
#'        e.g., generated from calls to \code{\link{randompar}}.
#' @param nsample Sample size
#' @return Data frame of response and predictor variables.
#' @export
randomsamp <- function(predictors, response, prob=NULL, par, nsample=100){
  
  Ly <- length(response)
  if(is.null(prob)) prob <- rep(1,Ly)
  if(length(prob)!=Ly) stop('prob length does not match nsample')
  y <- sample(response, size=nsample, replace=TRUE, prob=prob)
  dat <- NULL
  for(iy in seq_len(Ly)){
    ny <- sum(y==response[iy])
    xi <- sample_xi(nsample=ny, predictors=predictors,
                          h=par[[iy]]$h, J=par[[iy]]$J)
    yx <- cbind(data.frame(y=rep(response[iy],ny)),xi)
    dat <- rbind(dat, yx)
  }

  return(dat)
}