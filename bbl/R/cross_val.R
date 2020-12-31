#' Cross-Validation of BB Learning
#' 
#' Run multiple fittings of \code{bbl} model with training/validation
#' division of data
#' 
#' The \code{data} slot of \code{object} is split into training and validation 
#' subsets of (\code{nfold}-1):1 ratio. The model is trained with the
#' former and validated on the latter. Individual division/fold results are 
#' combined into validation result for all instances in the data set and
#' prediction score is evaluated using the known response group
#' identity.
#' 
#' @param formula Formula for model. Note that intercept has no effect.
#' @param data Data frame of data. Column names must match \code{formula}.
#' @param weights Frequency vector of how many times each row of \code{data} must
#'        be repeated. If \code{NULL}, defaults to vector of 1s. 
#'        Fractional values are not supported.
#' @param novarOk Proceed even when there are predictors with only one factor level.
#' @param lambda Vector of L2 penalizer values for \code{method = 'pseudo'}. Inferences
#'        will be repeated for each value. Restricited to non-negative values.
#' @param lambdah L2 penalizer in \code{method = 'pseudo'} applied to
#'        parameter \code{h}. In contrast to \code{lambda}, 
#'        only a single value is allowed.
#' @param eps Vector of regularization parameters, \eqn{\epsilon\in[0,1]}, 
#'        for \code{method = 'mf'}. Inference will be repeated
#'        for each value.
#' @param nfold Number of folds for training/validation split.
#' @param method \code{c('pseudo','mf')} for pseudo-likelihood maximization or
#'        mean field.
#' @param use.auc Use AUC as the measure of prediction accuracy. Only works
#'        if response groups are binary. If \code{FALSE}, mean prediction group
#'        accuracy will be used as score.
#' @param verbose Verbosity level. Downgraded when relayed into \code{\link{bbl}}.
#' @param progress.bar Display progress bar in \code{\link{predict}}.
#' @param storeOpt Store the optimal fitted object of class \code{\link{bbl}}.
#' @param ... Other parameters to \code{\link{mlestimate}}.
#' @return Object of class \code{cv.bbl} extending \code{\link{bbl}}, a list
#'         with extra components
#'         \item{regstar}{Value of regularization parameter, \code{lambda}
#'         and \code{eps} for \code{method='pseudo'} and \code{method='mf'},
#'         respectively, at which the accuracy score is maximized}
#'         \item{maxscore}{Value of maximum accuracy score}
#'         \item{cvframe}{Data frame of regularization parameters and scores scanned.
#'             If \code{use.auc=TRUE}, also contains 95% c.i.}
#'         The components of \code{\link{bbl}} store the optimal model trained
#'         if \code{storeOpt=TRUE}.
#' @examples
#' set.seed(513)
#' m <- 5
#' n <- 100
#' predictors <- list()
#' for(i in 1:m) predictors[[i]] <- c('a','c','g','t')
#' names(predictors) <- paste0('v',1:m)
#' par <- list(randompar(predictors), randompar(predictors, h0=0.1, J0=0.1))
#' dat <- randomsamp(predictors, response=c('ctrl','case'), par=par, nsample=n)
#' cv <- crossVal(y ~ .^2, data=dat, method='mf', eps=seq(0.1,0.9,0.1))
#' cv
#' @export
crossVal <- function(formula, data, weights=NULL, novarOk=FALSE, 
                     lambda=1e-5, lambdah=0, eps=0.9, nfold=5, 
                     method='pseudo', use.auc=TRUE, 
                     verbose=1, progress.bar=FALSE, storeOpt=TRUE, ...){
  
  
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
# data <- as.data.frame(lapply(data, function(x) if(is.numeric(x)) factor(x)))
  term <- terms(formula, data=data)
  idy <- attributes(term)$response
  vars <- as.character(attributes(term)$variables)
  colnames(data) <- fixnames(colnames(data))
  resp <- vars[[idy+1]]
  if(!resp %in% colnames(data)) stop('Response var not in data')
# vars <- vars[!(vars %in% c('list',resp))]
  vars <- vars[-1]   # 'list' can be in colnames 
  vars <- vars[!vars==resp]
# xlevels <- .getXlevels(term, m=data)
  xlevels <- getxlevels(vars, data=data)
  formula <- formula(term)
  for(i in seq_along(xlevels)){
    if(length(xlevels[[i]])==1){
      mesg <- paste0('Predictor ',names(xlevels)[i],' has one level')
      if(novarOk)
        warning(mesg)
      else
        stop(mesg)
    }
  }
  m <- length(xlevels)

  y <- data[,resp]
  yx <- data[,c(resp,names(xlevels))]
  groups <- levels(factor(y))
  Ly <- length(groups)
  
  if(Ly==1) warning('Only one response group in data')
  if(Ly!=2) use.auc <- FALSE
  
  if(!is.null(weights)){
    if(!all.equal(weights, as.integer(weights))) stop('Non-integer weights')
    if(length(weights)!=NROW(data))
      stop('Length of weights does not match data')
    yx <- freq2raw(data=yx, freq=weights)
    y <- yx[,resp]
  }
  
  if(method=='pseudo'){ 
    reglist <- lambda
    if(length(lambdah)>1) 
      stop('Only a single value of lambdah allowed')
  }
  else if(method=='mf') reglist <- eps

  else stop('Unknown method')
  
  res <- NULL    # data frame of cv result
  bbopt <- NULL  # optimally trained bbl object
  maxscore <- -Inf
  regstar <- 0
  for(reg in reglist){
    if(is.null(lambdah)) regh <- reg
    else regh <- lambdah
    if(verbose>0){ 
      if(method=='pseudo'){
        if(regh>0)
          cat('Cross validation under lambda = (',reg,
            ',',regh,')\n',sep='')
        else
          cat('Cross validation under lambda = ',reg,'\n',sep='')
      }
      else cat('Cross validation under epsilon = ',reg,'\n',sep='')
    }
    pred <- NULL
    for(k in seq_len(nfold)){
      if(verbose>0) cat(' Fold no. ',k,'...\n',sep='')
      itrain <- ival <- NULL
      for(iy in seq_len(Ly)){
        idy <- which(y==groups[iy])
        ns <- length(idy)
        nval <- max(1,floor(ns/nfold))
        imax <- k*nval
        if(imax > ns | k==nfold) imax <- ns
        if(imax < 1) imax <- 1
        iyval <- idy[seq((k-1)*nval+1, imax)]
        iytrain <- idy[!idy %in% iyval]
        ival <- c(ival, iyval)
        itrain <- c(itrain, iytrain)
      }
      if(sum(is.na(ival)>0) | sum(is.na(itrain)>0)) 
        stop('Error in train/validation division')
      if(length(ival)==0 | length(itrain)==0)
        stop('Sample size too small for nfold requested')
      dval <- yx[ival, ,drop=FALSE]
      dtrain <- yx[itrain, ,drop=FALSE]
      if(method=='pseudo')
        obtrain <- bbl(formula, data=dtrain, xlevels=xlevels, 
                       method=method, novarOk=TRUE, 
                       testNull=FALSE, lambda=reg, 
                       verbose=verbose-1, lambdah=regh, ...)
      else
        obtrain <- bbl(formula, data=dtrain, xlevels=xlevels,
                       method=method, novarOk=TRUE, 
                       testNull=FALSE, eps=reg, 
                       verbose=verbose-1, lambdah=lambdah, ...)

      pr <- predict(object=obtrain, newdata=dval, logit=!use.auc,
                    progress.bar=progress.bar, verbose=verbose-1)
      yxv <- data.frame(ytrue=factor(yx[ival,1],levels=groups))
      pred <- rbind(pred, cbind(yxv, pr))
    }
    if(use.auc){
      roc <- pROC::roc(response=pred[,1], levels=groups, 
                 predictor=pred[,3], direction='<', ci=TRUE)
      score <- roc$auc
      ci <- roc$ci
      if(verbose>0) cat(' AUC = ',score,'\n',sep='')
    }
    else{
      score <- mean(pred[,1]==pred$yhat)
      if(verbose>0) cat(' Prediction score = ',score,'\n',sep='')
      ci <- NULL
    }
    if(method=='pseudo'){
      if(use.auc) rx <- data.frame(lambda=reg, AUC=score, ci1=ci[1],
                                   ci2=ci[3])
      else rx <- data.frame(lambda=reg, score=score)
    } else{
      if(use.auc) rx <- data.frame(epsilon=reg, AUC=score, ci1=ci[1],
                                   ci2=ci[3])
      else rx <- data.frame(epsilon=reg, score=score)
    }
    if(score > maxscore){
      maxscore <- score
      regstar <- reg
      if(storeOpt) bbopt <- obtrain
    }
    res <- rbind(res, rx)
  }
  
  cv <- c(bbopt, list(regstar=regstar, maxscore=maxscore, cvframe=res))
  class(cv) <- 'cv.bbl'
  
  return(cv)
}

getxlevels <- function(vars, data){
  
  m <- length(vars)
  xlevels <- vector('list',m)
  names(xlevels) <- vars
  for(i in seq_along(vars)){
    if(!vars[i] %in% colnames(data)) 
      stop(paste0("'",vars[i],"' not in data"))
    v <- data[,vars[i]]
    if(is.factor(v) | is.character(v)){
      lv <- levels(factor(v))
      xlevels[[i]] <- lv[order(lv)]
    }
    else{
      v <- unique(v)
      xlevels[[i]] <- v[order(v)]
    }
  }
  return(xlevels)
}

# fix numeric/special char column names
fixnames <- function(cdat){

  for(i in seq_along(cdat)){
    x <- cdat[i]
    if(substring(x,1,1)=='`' & substring(x,nchar(x),nchar(x))=='`') next()  
    if(make.names(x)!=x)
      cdat[i] <- paste0('`',x,'`')
  }
  
  return(cdat)
}  