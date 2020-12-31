#' Boltzmann Bayes Learning Inference
#' 
#' Main driver for bbl inference
#' 
#' Formula argument and data are used to tabulate xlevels unless explicitly
#' given as list. Data are expected to be factors or integers. This function
#' is a driver interepreting formula and calls \code{bbi.fit}. Will stop with
#' error if any predictor has only one level unless \code{novarOk='TRUE'}.
#' Use \code{\link{removeConst}} to remove the non-varying predictors before 
#' calling if this happens.
#' 
#' @param formula Formula for modeling
#' @param data Data for fitting
#' @param weights Vector of weights for each instance in data. Restricted to 
#'        non-negative integer frequencies, recoding the number of times 
#'        each row of data must be repeated. If \code{NULL},
#'        assumed to be all 1. Fractional weights are not supported.
#' @param xlevels List of factor levels for predictors. If \code{NULL},
#'        will be inferred from data with factor levels ordered alphanumerically.
#' @param verbose Output verbosity level. Will be send to down-stream function
#'        calls with one level lower.
#' @param method BB inference algorithm; pseudo-likelihood inference (\code{'pseudo'})
#'        or mean field (\code{'mf'}).
#' @param novarOk If \code{TRUE}, will proceed with predictors having only one
#'        level.
#' @param testNull Repeat the inference for the `pooled' sample; i.e., under the
#'        null hypothesis of all rows in data belonging to a single group.
#' @param prior.count Prior count for computing single predictor and pairwise
#'        frequencies
#' @param ... Other parameters to \code{\link{mlestimate}}.
#' @return
#' A list of class \code{bbl} with the following elements:
#'   \item{coefficients}{List of inferred coefficients with elements
#'   \code{h}, \code{J}, \code{h0}, and \code{J0}. The bias  
#'   parameter \code{h} is a list of length equal to no. of 
#'   response groups, each of which is a list of the same struture as 
#'   \code{xlevels}: length equal to no. of predictors, containing vectors of 
#'   length equal to each predictor factor levels:
#'   \eqn{h_i^{(y)}(x)} represented by \code{h[[y]][[i]][x]}.
#'   The interaction parameter \code{J} is a list of lists of dimension
#'   \eqn{m \times m}, where \eqn{m} is the number of predictors. Each
#'   element is a matrix of dimension \eqn{L_i \times L_j}, where \eqn{L_i}
#'   and \eqn{L_j} are numbers of factor levels in predictor \code{i} and
#'   \code{j}: \eqn{J_{ij}^{(y)}(x_1,x_2)} represented by 
#'   \code{J[[y]][[i]][[j]][x1,x2]}. All elements of lists are named.
#'   The pooled parameters \code{h0} and \code{J0}, if computed,
#'   are of one less dimension, omitting response group argument.}
#'   \item{xlevels}{List of vectors containing predictor levels.}
#'   \item{terms}{The \code{terms} of \code{formula} input.}
#'   \item{groups}{Vector of response groups.}
#'   \item{groupname}{Name of the response variable.}
#'   \item{qJ}{Matrix of logicals whose elements record whether
#'     \code{formula} includes interaction between the two predictors.}
#'   \item{model}{Model data frame derived from \code{formula} and \code{data}.}
#'   \item{lkh}{Log likelihood.}
#'   \item{lz}{Vector log partition function. Used in \code{\link{predict}}.}
#'   \item{weights}{Vector of integral weights (frequencies).}
#'   \item{call}{Function call.}
#'   \item{df}{Degrees of freedom.} 
#' @examples
#' titanic <- as.data.frame(Titanic)
#' b <- bbl(Survived ~ .^2, data=titanic[,1:4], weights=titanic$Freq)
#' b
#' @import Rcpp
#' @import stats
#' @import graphics
#' @importFrom grDevices colorRampPalette gray.colors
#' @useDynLib bbl
#' @export
bbl <- function(formula, data, weights=NULL, xlevels=NULL, verbose=1, 
                method='pseudo', novarOk=FALSE, testNull=TRUE, 
                prior.count=1, ...){

  cl <- match.call()
  if(missing(data)) 
    stop('data argument required')
  if(!is.null(weights)){
    if(length(weights)!=NROW(data))
      stop('Length of weights does not match data')
    zero <- weights==0
    data <- data[!zero,]
    weights <- weights[!zero]   # remove rows with zero weights
  }
  
  term <- stats::terms(formula, data=data)
  idy <- attributes(term)$response
  vars <- as.character(attributes(term)$variables)
  colnames(data) <- fixnames(colnames(data))
  resp <- vars[[idy+1]]
# vars <- vars[!(vars %in% c('list',resp))]
  vars <- vars[-1]  # 'list' can be in colnames
  vars <- vars[!vars==resp]
# xlevels <- .getXlevels(term, m=data)
  if(is.null(xlevels))
    xlevels <- getxlevels(vars, data=data)
  
  formula <- formula(term)
  
  if(!novarOk){
    for(i in seq_along(xlevels)){
      if(length(xlevels[[i]])==1)
        stop(paste0('Predictor ',names(xlevels)[i],' has one level'))
    }
  }
  m <- length(xlevels)
  idy <- attributes(term)$response
  resp <- as.character(attributes(term)$variables[[idy+1]])
  y <- data[,resp]
  x <- data[,names(xlevels)]
  
  Ly <- length(levels(factor(y)))
  if(Ly==1) warning('Only one response group in data')
  
  label <- attr(term,'term.labels')
  ilabel <- label[vapply(label,FUN=function(x){grepl(':',x)}, logical(1))]
# ilabel <- gsub('`','',ilabel)
  ijlabel <- strsplit(ilabel,split=':')
  qJ <- matrix(FALSE, nrow=m, ncol=m)
  rownames(qJ) <- colnames(qJ) <- names(xlevels)
  for(k in seq_along(ijlabel))
    qJ[ijlabel[[k]][1],ijlabel[[k]][2]] <- TRUE
  qJ <- qJ | t(qJ)    # TRUE for all interacting pairs of predictors
  naive <- sum(qJ)==0
  groups <- levels(factor(y))
  
  bb <- list()
  class(bb) <- 'bbl'
  
  if(naive & method=='mf'){
    b <- naivemf(xlevels=xlevels, y=y, weights=weights, data=data, 
                 prior.count=prior.count)
    lz <- rep(0, Ly)
    for(iy in seq_len(Ly)){
      for(i in seq_len(m)) 
        lz[[iy]] <- lz[[iy]] + log(1+sum(exp(b$h[[iy]][[i]])))
      lz[[iy]] <- lz[[iy]]/m
    }
    bb$coefficients <- list(h=b$h, h0=b$h0)
    bb$lz <- lz
  } else{
    b <- bbl.fit(x=x, y=y, qJ=qJ, weights=weights, xlevels=xlevels,
               verbose=verbose-1, method=method, prior.count=prior.count,
               ...) # alternative
    if(testNull)
      b0 <- bbl.fit(x=x, y=rep('pooled',length(y)), qJ=qJ, weights=weights, 
                xlevels=xlevels, verbose=verbose-1, method=method,
                prior.count=prior.count, ...) # null
    else b0 <- NULL
    bb$coefficients <- list(h=b$h, J=b$J, h0=b0$h[[1]], J0=b0$J[[1]])
    bb$lz <- b$lz
  }
  bb$xlevels <- xlevels
  bb$terms <- term
  bb$groups <- levels(factor(y))
  bb$groupname <- resp
  bb$qJ <- qJ
  bb$model <- data[,c(resp, names(xlevels))]
  bb$lkh <- b$lkh
  bb$weights <- weights
  bb$call <- cl
  if(naive) bb$method <- 'mf'
  else bb$method <- method
  df <- 0
  for(i in seq_len(m)){
    ni <- length(xlevels[[i]])-1
    df <- df + ni
    if(naive) next()
    if(i==m) next()
    for(j in seq(i+1,m))
      if(qJ[i,j]) df <- df + ni*(length(xlevels[[j]])-1)
  }
  bb$df <- df
  
  return(bb)
}

#' bbl Inference with model matrix
#' 
#' Performs bbl inference using response vector and predictor matrix
#' 
#' This function would normally be called by \code{\link{bbl}} rather than
#' directly. Expects the predictor data \code{x} and response vector \code{y}
#' instead of formula input to \code{\link{bbl}}.
#' 
#' @param x Data frame of factors with each predictor in columns.
#' @param y Vector of response variables.
#' @param qJ Matrix of logicals indicating which predictor combinations
#'           are interacting. 
#' @param weights Vector of non-negative integer frequencies, recoding 
#'        the number of times each row of data must be repeated. 
#'        If \code{NULL}, assumed to be all 1. Fractional weights are not 
#'        supported.
#' @param xlevels List of factor levels for predictors. If \code{NULL},
#'        will be inferred from data with factor levels ordered alphanumerically.
#' @param verbose Verbosity level of output. Will be propagated to 
#'        \code{\link{mlestimate}} with one level down.
#' @param method \code{c('pseudo','mf')}; inference method.
#' @param prior.count Prior count for computing single predictor and pairwise
#'        frequencies
#' @param ... Other arguments to \code{\link{mlestimate}}.
#' @return List of named components \code{h}, \code{J}, \code{lkh}, and
#'         \code{lz}; see \code{\link{bbl}} for information regarding these
#'         components.
#' @examples
#' titanic <- as.data.frame(Titanic)
#' freq <- titanic$Freq
#' x <- titanic[,1:3]
#' y <- titanic$Survived
#' b <- bbl.fit(x=x,y=y, weights=freq)
#' b
#' @export
bbl.fit <- function(x, y, qJ=NULL, weights=NULL, xlevels=NULL, verbose=1, 
                    method='pseudo', prior.count=1, ...){
  
  le <- levels(factor(y))
  Ly <- length(le)
  h <- J <- list()
  m <- NCOL(x)
  vars <- colnames(x)
  if(is.null(xlevels))
    xlevels <- getxlevels(colnames(x), data=x)
  if(is.null(qJ)){
    qJ <- matrix(TRUE, nrow=m, ncol=m)
    diag(qJ) <- FALSE
    rownames(qJ) <- colnames(qJ) <- vars
  }
  if(!all.equal(dim(qJ), c(m,m)))
    stop('Incorrect dimension of qJ')
                       
  lkh <- 0
  lz <- rep(0, Ly)
  for(iy in seq_len(Ly)){
    ny <- sum(y==le[iy])
    if(ny==0) 
      stop(paste0('No instance of "',le[iy],'" in training data'))
    da <- x[y==le[iy],]
    if(!is.null(weights))
      frq <- weights[y==le[iy]]
    else frq <- rep(1, sum(y==le[iy]))
    if(verbose > 1) cat(' Inference for y = ',le[iy],':\n',sep='')
#   xda <- data.matrix(da) - 1
    xda <- matrix(0, nrow=NROW(da), ncol=NCOL(da))
    L <- NULL
    for(i in seq_len(NCOL(da))){
      xda[,i] <- match(da[,i],xlevels[[i]])-1
      L[i] <- length(xlevels[[i]])
    }
    b <- mlestimate(xi=xda, weights=frq, qJ=qJ, verbose=verbose-1, 
                    method=method, L=L, prior.count=prior.count, ...)
    names(b$h) <- names(b$J) <- names(xlevels)
    for(i in seq_len(m)){
      ni <- xlevels[[i]][-1][seq_along(b$h[[i]])]
      names(b$h[[i]]) <- ni
      names(b$J[[i]]) <- names(xlevels)
      for(j in seq_len(m)){
        if(is.null(b$J[[i]][[j]])) next()
        if(NROW(b$J[[i]][[j]])>0)
          rownames(b$J[[i]][[j]]) <- ni[seq_len(NROW(b$J[[i]][[j]]))]
        if(NCOL(b$J[[i]][[j]])>0)
          colnames(b$J[[i]][[j]]) <- 
              xlevels[[j]][-1][seq_len(NCOL(b$J[[i]][[j]]))]
      }
    }
    h[[iy]] <- b$h
    J[[iy]] <- b$J
    lkh <- lkh + b$lkh*ny
    lz[iy] <- b$lz
  }
  names(lz) <- names(h) <- names(J) <- le

  return(list(h=h, J=J, lkh=lkh, lz=lz))
}
