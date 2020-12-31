binomial.fit.robmixglm <- function(x,y,offset,gh,notrials,EMTol,  calcHessian=TRUE,  cores, verbose,starting.values=NULL) {
  
  if (!is.null(starting.values)) notrials <- 1
  
  logit <- function(x){
    1.0/(1.0+exp(-x))
  }
  
  llrandbinom <- function(y, lp, tau2, gh){
    
    dointegrate <- function(oney,onelp,tau2) {
      
      thefun <- function(nu){
        return(dbinom(oney[1],oney[1]+oney[2],logit(onelp+nu*sqrt(tau2))))
      }
      theint <- sum(thefun(gh[,1])*gh[,2])
      if (is.nan(theint)) theint <- NA
      return(log(theint))
    }
    
    calconell <- function(x) {
      oney <- x[1:2]
      onelp <- x[3]
      if (tau2<0.0) return(NA)
      if (tau2==0) thell <- dbinom(oney[1],oney[1]+oney[2],logit(onelp),log=TRUE) else
        thell <- dointegrate(oney=oney,onelp=onelp,tau2=tau2)
      return(thell)
    }
    if (is.na(tau2)) stop("tau2 is NA")
    ll <- apply(cbind(y,lp),1,calconell)
    return(ll)
  }
  
  fitonemlreg <- function(y,outliers,x=NULL,offset,fixed) {
    
    rlreg <- function(xcoef,lpoutlier,tau2,x,y,offset,prop) {
      poutlier <- 1.0/(1+exp(-lpoutlier))
      
      lp <- as.vector(x %*% xcoef)+offset
      
      lp1 <- dbinom(y[,1], y[,1]+y[,2], logit(lp), log = TRUE)
      lp2 <- llrandbinomcpp(y, lp, tau2, gh)
      
      if (!missing(prop)) {
        ll <- prop*cbind(lp1,lp2)
        negll <- -sum(apply(ll,1,sum))
      } else {
        l <- exp(cbind(lp1+log(1-poutlier),lp2+log(poutlier)))
        negll <- -sum(log(apply(l,1,sum)))
      }
      if (is.nan(negll)) negll <- 1e180
      if (!is.finite(negll)) negll <- 1e180
      return(negll)
    }
    
    optimrlreg <- function(p,lpoutlier,x,y,offset,prop,fixed) {
      p[names(fixed)] <- fixed
      return(rlreg(matrix(p[1:(length(p)-1)],ncol=1),lpoutlier,p[length(p)],x,y,offset,prop))
    }
    
    tryCatch({
      if (is.null(starting.values)) {
        robust.binomial.prefit <- glm(cbind(y[,1],y[,2])~x[,colnames(x)!="(Intercept)"],family=binomial(), offset=offset,
                                      subset=(outliers!=1))
        
        
        # get the starting values
        prefit.coef <- coef(robust.binomial.prefit)
        # assume 20% outliers as a starting point
        currlpoutlier <- log(0.2/(1-0.2))
        currxcoef <- matrix(prefit.coef[1:(length(prefit.coef))],ncol=1)
        currxcoef <- ifelse(is.na(currxcoef),0,currxcoef)
        currtau2 <- min(rgamma(1,2,1),2)
      } else
      {
        currxcoef <- matrix(starting.values[1:(length(starting.values)-2)],ncol=1)
        currlpoutlier <- starting.values[length(starting.values)-1]
        currtau2 <- starting.values[length(starting.values)]
      }
      
      currll <- -1.0e100
      nem <- 0
      
      repeat {
        nem <- nem+1
        # expectation step
        currpoutlier <- 1.0/(1+exp(-currlpoutlier))
        
        lp <- as.vector(x %*% currxcoef)+offset
        
        ll1 <- dbinom(y[,1], y[,1]+y[,2], logit(lp), log = TRUE)+log(1-currpoutlier) 
        ll2 <- llrandbinomcpp(y, lp, currtau2, gh)+log(currpoutlier)
        
        ll <- cbind(ll1,ll2)
        prop <- t(apply(ll,1,function(x) {
          x <- x-max(x)
          x <- ifelse(x==-Inf,-1e100,x)
          return(exp(x)/sum(exp(x)))
        }))
        # calculate outlier proportion
        poutlier <- sum(prop[,2])/dim(prop)[1]
        currlpoutlier <- log(poutlier/(1-poutlier))
        
        if (is.na(poutlier)) stop("poutlier is NA")
        
        startvals <- c(currxcoef,currtau2)
        names(startvals) <- c(dimnames(x)[[2]],"tau2")
        
        results.nlm <- suppressWarnings(nlminb(startvals,optimrlreg,
                                               lower=c(rep(-Inf,length(currxcoef)),0),
                                               #control=list(trace=1,iter.max=10),
                                               control=list(iter.max=5),
                                               lpoutlier=currlpoutlier,
                                               prop=prop,
                                               y=y,x=x,offset=offset,
                                               fixed=fixed))
        currxcoef <- matrix(as.numeric(results.nlm$par)[1:(length(results.nlm$par)-1)],ncol=1)
        currtau2 <- as.numeric(results.nlm$par)[length(results.nlm$par)]
        
        lastll <- currll
        currll <- -rlreg(currxcoef,currlpoutlier,currtau2,x,y,offset)
        if (verbose) print(sprintf("Likelihood at end of EM step %.4f", currll))
        if (abs((lastll-currll)/currll)<EMTol) break()
        if (nem >100) break()
      }
      return(list(ll=currll,start.val=c(currxcoef,currlpoutlier,currtau2)))    
    },
    error=function(e) return(list(ll=NA))
    )
  }
  
  ll.robustbinomial <- function(p){
    
    xcoef <- p[1:(length(p)-2)]
    lpoutlier <- p[length(p)-1]
    tau2 <- p[length(p)]
    
    poutlier <- exp(lpoutlier)/(1+exp(lpoutlier))
    
    lp <- as.vector(x %*% xcoef)+offset
    
    ll1 <- dbinom(y[,1], y[,1]+y[,2], logit(lp), log = TRUE)+log(1-poutlier)
    
    ll2 <- llrandbinomcpp(y, lp, tau2, gh)+log(poutlier)
    
    l <- exp(cbind(ll1,ll2))
    negll <- -sum(log(apply(l,1,sum)))
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
  
  thenames <- c(dimnames(x)[[2]],"lpoutlier","tau2")
  
  names(start.val) <- thenames
  
  parnames(ll.robustbinomial) <- names(start.val)
  
  lower.val <- c(rep(-Inf,length(start.val)-2),-Inf,0)
  
  names(lower.val) <- names(start.val)
 
   if(verbose) thecontrol <- list(eval.max=1000,iter.max=1000,eval.max=2000,trace=5)
  else thecontrol <- list(eval.max=1000,iter.max=1000,eval.max=2000)
 
  robustbinomial.fit <- mle2(ll.robustbinomial,start=start.val,vecpar=TRUE,
                             optimizer="user",optimfun=myoptim,
                             data=list(y=y,x=x,offset=offset),
                             skip.hessian=TRUE,trace=verbose,
                             lower=lower.val,
                            control=thecontrol)
  if (calcHessian) {
    thecoef <- coef(robustbinomial.fit)
    ncoef <- length(thecoef)
    if (thecoef[ncoef]<0.0001) thecoef[ncoef] <- 0.0001
    robustbinomial.fit@details$hessian <- optimHess(thecoef,ll.robustbinomial,control=list(ndeps=c(rep(1.0e-5,length(thecoef)))))
    robustbinomial.fit@vcov <- ginv(robustbinomial.fit@details$hessian)
  }
  # if (any(is.nan(sqrt(diag(robustbinomial.fit@vcov)))))  warning("Error in calculating standard errors.")
  
  xcoef <- matrix(coef(robustbinomial.fit)[1:(length(coef(robustbinomial.fit))-2)],ncol=1)
  lpoutlier <- coef(robustbinomial.fit)[length(coef(robustbinomial.fit))-1]
  poutlier <- 1.0/(1+exp(-lpoutlier))
  tau2 <- coef(robustbinomial.fit)[length(coef(robustbinomial.fit))]
  
  lp <- as.vector(x %*% xcoef)+offset
  
  ll1 <- dbinom(y[,1], y[,1]+y[,2], logit(lp), log = TRUE)+log(1-poutlier) 
  ll2 <- llrandbinomcpp(y, lp, tau2, gh)+log(poutlier)
  
  ll <- cbind(ll1,ll2)
  prop <- t(apply(ll,1,function(x) {
    x <- x-max(x)
    x <- ifelse(x==-Inf,-1e100,x)
    return(exp(x)/sum(exp(x)))
  }))
  
  coef.names <- c(dimnames(x)[[2]],"Outlier p.","Tau-sq")
  return(list(fit=robustbinomial.fit,prop=prop,logLik=-robustbinomial.fit@min,np=length(coef.names),nobs=dim(x)[1],coef.names=coef.names))  
}

