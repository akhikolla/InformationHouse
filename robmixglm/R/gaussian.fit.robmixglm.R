gaussian.fit.robmixglm <- function(x,y,offset,gh,notrials,EMTol,  calcHessian=TRUE, cores, verbose,   starting.values=NULL) {
  
  if (!is.null(starting.values)) notrials <- 1
  
  fitonemlreg <- function(y,outliers,x=NULL,offset=NULL,fixed) {
    
    rlreg <- function(xcoef,lpoutlier,tau2, sigma2,x,y,offset,prop) {
      poutlier <- 1.0/(1+exp(-lpoutlier))
      
      lp <- as.vector(x %*% xcoef)+offset
      
      lp1 <- dnorm(y, mean=lp, sd=sqrt(sigma2), log = TRUE)
      lp2 <- dnorm(y, mean=lp, sd=sqrt(tau2+sigma2), log=TRUE)
      
      if (!missing(prop)) {
        ll <- prop*cbind(lp1,lp2)
        negll <- -sum(apply(ll,1,sum))
      } else {
        l <- exp(cbind(lp1+log(1-poutlier),lp2+log(poutlier)))
        negll <- -sum(log(apply(l,1,sum)))
      }
      if (is.nan(negll)) negll <- NA
      if (!is.finite(negll)) negll <- NA
      return(negll)
    }
    
    optimrlreg <- function(p,lpoutlier,x,y,offset,prop,fixed) {
      p[names(fixed)] <- fixed
      return(rlreg(matrix(p[1:(length(p)-2)],ncol=1),lpoutlier,p[length(p)-1],p[length(p)],x,y,offset,prop))
    }
          
    tryCatch({
      if (is.null(starting.values)) {
      robust.gaussian.prefit <- lm(y~x[,colnames(x)!="(Intercept)"],offset=offset,subset=(outliers!=1))
      prefit.coef <- coef(robust.gaussian.prefit)
      # assume 20% outliers as a starting point
      currlpoutlier <- log(0.2/(1-0.2))
      currxcoef <- matrix(prefit.coef[1:(length(prefit.coef))],ncol=1)
      currxcoef <- ifelse(is.na(currxcoef),0,currxcoef)
      currtau2 <- min(rgamma(1,2,1),5)*summary(robust.gaussian.prefit)$sigma^2
      currsigma2 <- summary(robust.gaussian.prefit)$sigma^2
    } else {
      currxcoef <- matrix(starting.values[1:(length(starting.values)-3)],ncol=1)
      currlpoutlier <- starting.values[length(starting.values)-2]
      currtau2 <- starting.values[length(starting.values)-1]
      currsigma2 <- starting.values[length(starting.values)]
    }
    currll <- -1.0e100
    nem <- 0
    
    repeat {
      nem <- nem+1
      # expectation step
      currpoutlier <- 1.0/(1+exp(-currlpoutlier))
      
      lp <- as.vector(x %*% currxcoef)+offset
      
      ll1 <- dnorm(y, mean=lp, sd=sqrt(currsigma2), log = TRUE)+log(1-currpoutlier) 
      ll2 <- dnorm(y, mean=lp, sd=sqrt(currsigma2+currtau2),log=TRUE)+log(currpoutlier)
      
      ll <- cbind(ll1,ll2)
      prop <- t(apply(ll,1,function(x) {
        x <- x-max(x)
        x <- ifelse(x==-Inf,-1e100,x)
        return(exp(x)/sum(exp(x)))
      }))
      # calculate outlier proportion
      poutlier <- sum(prop[,2])/dim(prop)[1]
      currlpoutlier <- log(poutlier/(1-poutlier))
      
      if (is.na(poutlier)) stop()
      
      startvals <- c(currxcoef,currtau2,currsigma2)
      names(startvals) <- c(dimnames(x)[[2]],"tau2","sigma2")
      
      results.nlm <- suppressWarnings(nlminb(startvals,optimrlreg,
                                             lower=c(rep(-Inf,length(currxcoef)),rep(0,2)),
                                             #control=list(trace=1,iter.max=10),
                                             control=list(trace=0,iter.max=5),
                                             lpoutlier=currlpoutlier,
                                             prop=prop,
                                             y=y,x=x,offset=offset,
                                             fixed=fixed))
   
      currxcoef <- matrix(as.numeric(results.nlm$par)[1:(length(results.nlm$par)-2)],ncol=1)
      currtau2 <- as.numeric(results.nlm$par)[length(results.nlm$par)-1]
      currsigma2 <- as.numeric(results.nlm$par)[length(results.nlm$par)]
      
      lastll <- currll
      currll <- -rlreg(currxcoef,currlpoutlier,currtau2,currsigma2,x,y,offset)
      if (verbose) print(sprintf("Likelihood at end of EM step %.4f", currll))
       if (abs((lastll-currll)/currll)<EMTol) break()
      if (nem >100) break()
    }
    return(list(ll=currll,start.val=c(currxcoef,currlpoutlier,currtau2,currsigma2)))    
    },
    error=function(e) return(list(ll=NA))
    )
  }
  
  ll.robustgaussian <- function(p){
    
    xcoef <- p[1:(length(p)-3)]
    lpoutlier <- p[length(p)-2]
    tau2 <- p[length(p)-1]
    sigma2 <- p[length(p)]
    
    poutlier <- 1.0/(1+exp(-lpoutlier))
    
    lp <- as.vector(x %*% xcoef)+offset
    
    ll1 <- dnorm(y, mean=lp, sd=sqrt(sigma2), log = TRUE)+log(1-poutlier)
    
    ll2 <- dnorm(y, mean=lp, sd=sqrt(sigma2+tau2), log = TRUE)+log(poutlier)
    
    ll <- cbind(ll1,ll2)
    maxll <- apply(ll,1,max)
    #browser()
    negll <- -sum(maxll+log(apply(exp(ll-maxll),1,sum)))
    if (is.nan(negll)) negll <- NA
    if (!is.finite(negll)) negll <- NA
    return(negll)
  }
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) runif(1)
    seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  
  if (is.null(starting.values)) {
    if (cores > 1) {
      cl = parallel::makeCluster(cores, setup_strategy = "sequential")
      doParallel::registerDoParallel(cl)
      res = foreach(i = 1:notrials, 
                    .options.RNG=seed[1]) %dorng% {
                      noutliers <- max(1,round(dim(x)[1]*0.2))
                      outliers <- sample(c(rep(1,noutliers),rep(0,dim(x)[1]-noutliers)),dim(x)[1])
                      fitonemlreg(y,outliers,x,offset,fixed=NULL)}
      parallel::stopCluster(cl)
#      browser()
      maxll <- -Inf
      nfails <- 0
      for (i in 1:notrials) {
        if (verbose) cat(c(res[[i]]$ll,res[[i]]$start.val),"\n")
        if (is.na(res[[i]]$ll)) nfails <- nfails+1
        else {
          if (res[[i]]$ll>maxll) {
            maxll <- res[[i]]$ll
            start.val <- res[[i]]$start.val
          }
        }
      }
      if (nfails > 0) warning(sprintf("Failed to obtain starting values for %i starting sets", nfails))
    } else {
      maxll <- -Inf
      nfails <- 0
      for (i in 1:notrials) {
          noutliers <- max(1,round(dim(x)[1]*0.2))
          outliers <- sample(c(rep(1,noutliers),rep(0,dim(x)[1]-noutliers)),dim(x)[1])
          if (verbose) print(sprintf("Trial %i", i))
          thefit <- fitonemlreg(y,outliers,x,offset,fixed=NULL)
          if (verbose) print(sprintf("Likelihood for trial %.4f", thefit$ll))
          if (verbose) print(thefit$start.val)
          if (is.na(thefit$ll)) nfails <- nfails+1
          else {
            if (thefit$ll>maxll) {
              maxll <- thefit$ll
              start.val <- thefit$start.val
            }
          }
        }
      if (nfails > 0) warning(sprintf("Failed to obtain starting values for %i starting sets", nfails))
    }
  } else {
    start.val <- starting.values
  }
    
  thenames <- c(dimnames(x)[[2]],"lpoutlier","tau2","sigma2")
  
  names(start.val) <- thenames
  
  parnames(ll.robustgaussian) <- names(start.val)
  
  lower.val <- c(rep(-Inf,length(start.val)-3),-Inf,0,0)
  
  names(lower.val) <- names(start.val)
  
  
  if(verbose) thecontrol <- list(eval.max=1000,iter.max=1000,eval.max=2000,trace=5)
  else thecontrol <- list(eval.max=1000,iter.max=1000,eval.max=2000)
 
  robustgaussian.fit <- mle2(ll.robustgaussian,start=start.val,vecpar=TRUE,
                             optimizer="user",optimfun=myoptim,
                             data=list(y=y,x=x,offset=offset),
                             skip.hessian=TRUE,trace=verbose,
                             lower=lower.val,
                             control=thecontrol)
  
  if (calcHessian) {
    thecoef <- coef(robustgaussian.fit)
    ncoef <- length(thecoef)
    if (thecoef[ncoef-1]<1e-4) thecoef[ncoef-1] <- 1e-4
    robustgaussian.fit@details$hessian <- optimHess(thecoef,ll.robustgaussian,control=list(ndeps=c(rep(1.0e-5,length(thecoef)))))
    robustgaussian.fit@vcov <- ginv(robustgaussian.fit@details$hessian)
  }
  # robustgaussian.fit@vcov <- myvcov
  #if (any(is.nan(sqrt(diag(robustgaussian.fit@vcov))))) warning("Error in calculating standard errors.")
  
  xcoef <- matrix(coef(robustgaussian.fit)[1:(length(coef(robustgaussian.fit))-3)],ncol=1)
  lpoutlier <- coef(robustgaussian.fit)[length(coef(robustgaussian.fit))-2]
  poutlier <- 1.0/(1+exp(-lpoutlier))
  tau2 <- coef(robustgaussian.fit)[length(coef(robustgaussian.fit))-1]
  sigma2 <- coef(robustgaussian.fit)[length(coef(robustgaussian.fit))]
  
  
  lp <- as.vector(x %*% xcoef)+offset
  
  ll1 <- dnorm(y, mean=lp, sd=sqrt(sigma2),log = TRUE)+log(1-poutlier) 
  ll2 <- dnorm(y, mean=lp, sd=sqrt(sigma2+tau2),log = TRUE)+log(poutlier)
  
  ll <- cbind(ll1,ll2)
  prop <- t(apply(ll,1,function(x) {
    x <- x-max(x)
    x <- ifelse(x==-Inf,-1e100,x)
    return(exp(x)/sum(exp(x)))
  }))
  
  coef.names <- c(dimnames(x)[[2]],"Outlier p.","Tau-sq","Sigma-sq")
  return(list(fit=robustgaussian.fit,prop=prop,logLik=-robustgaussian.fit@min,np=length(coef.names),nobs=dim(x)[1],coef.names=coef.names))  
}

