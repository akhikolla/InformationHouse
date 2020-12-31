gaussian.outlierTest.robmixglm <- function(object, R, parallel, cores) {
  
  gaussian.rg <- function(data, mle) {
    
    # calculate linear predictor
    lp <- data$X %*% matrix(mle$coefficients,ncol=1) + data$offset
    sd <- sqrt(mle$deviance/mle$df.residual)
    data$Y <- rnorm(length(data$Y),mean=lp,sd=sd)
    return(data)
  }
  
  gaussian.fun <- function(thedata) {
    gaussian.mle2 <- glm.fit(thedata$X, thedata$Y, 
                             offset = thedata$offset, family = gaussian())
    thesd <- sqrt(gaussian.mle2$deviance/gaussian.mle2$df.residual)
    # starting values assume 20% outliers and tau^2 of 1
    starting.values <- c(coef(gaussian.mle2)[1:length(coef(gaussian.mle2))],log(0.2/(1-0.2)),1,thesd)
    if (fitno==1) mixfit <- suppressWarnings(gaussian.fit.robmixglm(thedata$X, thedata$Y, 
                                       offset = thedata$offset,gh = norm.gauss.hermite(object$quadpoints),
                                       notrials=object$notrials, EMTol = object$EMTol, calcHessian=FALSE,
                                       verbose=FALSE, starting.values=NULL, cores=1))
     else mixfit <- suppressWarnings(gaussian.fit.robmixglm(thedata$X, thedata$Y, 
                                       offset = thedata$offset,gh = norm.gauss.hermite(object$quadpoints),
                                       notrials=object$notrials, EMTol = object$EMTol, calcHessian=FALSE,
                                       verbose=FALSE, starting.values=starting.values, cores=1))
    lp <- thedata$X %*% matrix(gaussian.mle2$coefficients,ncol=1) + thedata$offset
    gaussian.loglik <- sum(dnorm(thedata$Y,mean=lp,sd=thesd,log=TRUE))
    sim.chisq <- 2*(mixfit$logLik-gaussian.loglik)
    fitno <<- fitno+1
    return(sim.chisq)
  }
  
  # fit glm to data
  gaussian.mle <- glm.fit(object$X, object$Y, 
                          offset = object$offset, family = gaussian())
  
  thedata <- list(X=object$X,Y=object$Y,offset=object$offset,mle=gaussian.mle)
  fitno <- 1
  gaussian.boot <- boot(thedata, gaussian.fun, R, sim="parametric",
                        ran.gen = gaussian.rg, mle=gaussian.mle,parallel=parallel,
                        ncpus = cores)
  return(gaussian.boot)
}

nbinom.outlierTest.robmixglm <- function(object, R, parallel, cores) {
  
  nbinom.rg <- function(data, mle) {
    
    # calculate linear predictor
    lp <- matrix(data$X,ncol=length(mle$par)-1) %*% matrix(mle$par[1:(length(mle$par)-1)],ncol=1) + data$offset
    theta <- mle$par[length(mle$par)]
    data$Y <- rnbinom(length(data$Y),mu=exp(lp),size=theta)
    return(data)
  }
  
  nbinom.fun <- function(data) {

    nbinom.mle2 <- fitnegbin(data$Y,object$X,data$offset)

    # starting values assume 20% outliers and tau^2 of 1
    starting.values <- c(nbinom.mle2$par[1:(length(nbinom.mle2$par)-1)],log(0.2/(1-0.2)),1,exp(nbinom.mle2$par[length(nbinom.mle2$par)]))
    if (fitno==1) mixfit <- suppressWarnings(nbinom.fit.robmixglm(data$X, data$Y, 
                                       offset = data$offset,gh = norm.gauss.hermite(object$quadpoints),
                                       notrials=object$notrials, EMTol = object$EMTol, calcHessian=FALSE,
                                       verbose=object$verbose, starting.values=NULL, cores=1))
     else mixfit <- suppressWarnings(nbinom.fit.robmixglm(data$X, data$Y, 
                                       offset = data$offset,gh = norm.gauss.hermite(object$quadpoints),
                                       notrials=object$notrials, EMTol = object$EMTol, calcHessian=FALSE,
                                       verbose=object$verbose, starting.values=starting.values, cores=1))
    

        lp <- data$X %*% matrix(nbinom.mle2$par[1:(length(nbinom.mle2$par)-1)],ncol=1) + data$offset
    
    nbinom.loglik <- sum(dnbinom(data$Y,mu=exp(lp),size=nbinom.mle2$par[length(nbinom.mle2$par)],log=TRUE))
    if (any(is.nan(dnbinom(thedata$Y,mu=exp(lp),size=nbinom.mle2$par[length(nbinom.mle2$par)],log=TRUE)))) browser()
    sim.chisq <- 2*(mixfit$logLik-nbinom.loglik)
    fitno <<- fitno+1
    return(sim.chisq)
  }

  nbinom.mle <- fitnegbin(object$Y,object$X,object$offset)
 
  thedata <- list(X=object$X,Y=object$Y,offset=object$offset,mle=nbinom.mle)
  fitno <- 1
  nbinom.boot <- boot(thedata, nbinom.fun, R, sim="parametric",
                        ran.gen = nbinom.rg, mle=nbinom.mle,parallel=parallel,
                      ncpus = cores)
  return(nbinom.boot)
}


gamma.outlierTest.robmixglm <- function(object, R, parallel, cores) {
  
  gamma.rg <- function(data, mle) {
    
    # calculate linear predictor
    lp <- data$X %*% matrix(mle$coefficients,ncol=1) + data$offset
    phi <- sum(mle$residuals^2)/mle$df.residual
     data$Y <- rgamma(length(data$Y),shape=1.0/phi,rate=1.0/(phi*exp(lp)))
    return(data)
  }
  
  gamma.fun <- function(thedata) {
    gamma.mle2 <- glm.fit(thedata$X, thedata$Y, 
                          offset = thedata$offset, family = Gamma(link="log"))
    # starting values assume 20% outliers and tau^2 of 1
    starting.values <- c(coef(gamma.mle2)[1:length(coef(gamma.mle2))],log(0.2/(1-0.2)),1,sum(gamma.mle2$residuals^2)/gamma.mle2$df.residual)
    if (fitno==1)  mixfit <- suppressWarnings(gamma.fit.robmixglm(thedata$X, thedata$Y, 
                                    offset = thedata$offset,gh = norm.gauss.hermite(object$quadpoints),
                                    notrials=object$notrials, EMTol = object$EMTol, calcHessian=FALSE,
                                    verbose=FALSE, starting.values=NULL, cores=1))
     else mixfit <- suppressWarnings(gamma.fit.robmixglm(thedata$X, thedata$Y, 
                                    offset = thedata$offset,gh = norm.gauss.hermite(object$quadpoints),
                                    notrials=object$notrials, EMTol = object$EMTol, calcHessian=FALSE,
                                    verbose=FALSE, starting.values=starting.values, cores=1))
    lp <- thedata$X %*% matrix(gamma.mle2$coefficients,ncol=1) + thedata$offset
    phi <- sum(gamma.mle2$residuals^2)/gamma.mle2$df.residual
    gamma.loglik <- sum(dgamma(thedata$Y,shape=1.0/phi,rate=1.0/(phi*exp(lp)),log=TRUE))
    sim.chisq <- 2*(mixfit$logLik-gamma.loglik)
     fitno <<- fitno+1
    return(sim.chisq)
  }
  
  # fit glm to data
  gamma.mle <- glm.fit(object$X, object$Y, 
                       offset = object$offset, family = Gamma(link="log"))
  
  
  thedata <- list(X=object$X,Y=object$Y,offset=object$offset,mle=gamma.mle)
  fitno <- 1
  gamma.boot <- boot(thedata, gamma.fun, R, sim="parametric",
                        ran.gen = gamma.rg, mle=gamma.mle, parallel=parallel,
                     ncpus = cores)
  return(gamma.boot)
}

truncpoisson.outlierTest.robmixglm <- function(object, R, parallel, cores) {
  
  truncpoisson.rg <- function(data, mle) {
    
    lp <- data$X %*% matrix(mle@coefficients,ncol=1) + data$offset
    data$Y <- actuar::rztpois(length(data$Y),exp(lp))
    return(data)
  }
  
  truncpoisson.fun <- function(data) {
      
    truncpoisson.mle2 <- vglm(data$Y~data$X[,colnames(data$X)!="(Intercept)"],
                                       family=pospoisson, offset=data$offset)
    
    # starting values assume 20% outliers and tau^2 of 1
    starting.values <- c(coef(truncpoisson.mle2)[1:length(coef(truncpoisson.mle2))],log(0.2/(1-0.2)),1)
    if (fitno==1) mixfit <- suppressWarnings(truncpoisson.fit.robmixglm(data$X, data$Y, 
                                           offset = thedata$offset,gh = norm.gauss.hermite(object$quadpoints),
                                           notrials=object$notrials, EMTol = object$EMTol, calcHessian=FALSE,
                                           verbose=object$verbose, starting.values=NULL, cores=1))
    else mixfit <- suppressWarnings(truncpoisson.fit.robmixglm(data$X, data$Y, 
                                           offset = thedata$offset,gh = norm.gauss.hermite(object$quadpoints),
                                           notrials=object$notrials, EMTol = object$EMTol, calcHessian=FALSE,
                                           verbose=object$verbose, starting.values=starting.values, cores=1))
    lp <- data$X %*% matrix(truncpoisson.mle2@coefficients,ncol=1) + data$offset
    truncpoisson.loglik <- sum(dztpois(data$Y,exp(lp),log=TRUE))
    sim.chisq <- 2*(mixfit$logLik-truncpoisson.loglik)
    fitno <<- fitno+1
    return(sim.chisq)
  }
  
  
  truncpoisson.mle2 <- vglm(object$Y~object$X[,colnames(object$X)!="(Intercept)"],family=pospoisson, offset=object$offset)
  
  fitno <- 1
  
  thedata <- list(X=object$X,Y=object$Y,offset=object$offset,mle=truncpoisson.mle2)
  
  truncpoisson.boot <- boot(thedata, truncpoisson.fun, R, sim="parametric",
                            ran.gen = truncpoisson.rg, mle=truncpoisson.mle2, parallel=parallel,
                            ncpus = cores)
  return(truncpoisson.boot)
}



poisson.outlierTest.robmixglm <- function(object, R, parallel, cores) {

    poisson.rg <- function(data, mle) {
      
      # calculate linear predictor
      lp <- data$X %*% matrix(mle$coefficients,ncol=1) + data$offset
      data$Y <- rpois(length(data$Y),exp(lp))
      return(data)
    }

    poisson.fun <- function(thedata) {
      poisson.mle2 <- glm.fit(thedata$X, thedata$Y, 
                              offset = thedata$offset, family = poisson())
      # starting values assume 20% outliers and tau^2 of 1
      starting.values <- c(coef(poisson.mle2)[1:length(coef(poisson.mle2))],log(0.2/(1-0.2)),1)
      if (fitno==1) mixfit <- suppressWarnings(poisson.fit.robmixglm(thedata$X, thedata$Y, 
                                        offset = thedata$offset,gh = norm.gauss.hermite(object$quadpoints),
                                        notrials=object$notrials, EMTol = object$EMTol, calcHessian=TRUE,
                                        verbose=object$verbose, starting.values=NULL, cores=1))
      else mixfit <- suppressWarnings(poisson.fit.robmixglm(thedata$X, thedata$Y, 
                                        offset = thedata$offset,gh = norm.gauss.hermite(object$quadpoints),
                                        notrials=object$notrials, EMTol = object$EMTol, calcHessian=FALSE,
                                        verbose=object$verbose, starting.values=starting.values, cores=1))
      lp <- thedata$X %*% matrix(poisson.mle2$coefficients,ncol=1) + thedata$offset
      poisson.loglik <- sum(dpois(thedata$Y,exp(lp),log=TRUE))
      sim.chisq <- 2*(mixfit$logLik-poisson.loglik)
    	fitno <<- fitno+1
      return(sim.chisq)
    }
    
    # fit glm to data
  poisson.mle <- glm.fit(object$X, object$Y, 
                         offset = object$offset, family = poisson())
  
 
  thedata <- list(X=object$X,Y=object$Y,offset=object$offset,mle=poisson.mle)
  fitno <- 1
  poisson.boot <- boot(thedata, poisson.fun, R, sim="parametric",
                       ran.gen = poisson.rg, mle=poisson.mle, parallel=parallel,
                       ncpus = cores)
  return(poisson.boot)
}

binomial.outlierTest.robmixglm <- function(object, R, parallel, cores) {
  
  binomial.rg <- function(data, mle) {
    # calculate linear predictor
    lp <- data$X %*% matrix(mle$coefficients,ncol=1) + data$offset
    n <- data$Y[,1]+data$Y[,2]
    data$Y[,1] <- rbinom(dim(data$Y)[1],data$Y[,1]+data$Y[,2],1.0/(1.0+exp(-lp)))
    data$Y[,2] <- n-data$Y[,1]
    return(data)
  }
  
  binomial.fun <- function(thedata) {
    binomial.mle2 <- glm.fit(thedata$X, thedata$Y,
                             offset = thedata$offset, family = binomial())
    starting.values <- c(coef(binomial.mle2),log(0.2/(1-0.2)),1)
    if (fitno==1) mixfit <- suppressWarnings(binomial.fit.robmixglm(thedata$X, thedata$Y,
                                       offset = thedata$offset,gh = norm.gauss.hermite(object$quadpoints),
                                       notrials=object$notrials, EMTol = object$EMTol, calcHessian=FALSE,
                                       verbose=FALSE, starting.values=NULL, cores=1))
    else mixfit <- suppressWarnings(binomial.fit.robmixglm(thedata$X, thedata$Y,
                                       offset = thedata$offset,gh = norm.gauss.hermite(object$quadpoints),
                                       notrials=object$notrials, EMTol = object$EMTol, calcHessian=FALSE,
                                       verbose=FALSE, starting.values=starting.values, cores=1))
    lp <- thedata$X %*% matrix(binomial.mle2$coefficients,ncol=1) + thedata$offset
    binomial.loglik <- sum(dbinom(thedata$Y[,1],thedata$Y[,1]+thedata$Y[,2],1.0/(1.0+exp(-lp)),log=TRUE))
     sim.chisq <- 2*(mixfit$logLik-binomial.loglik)
     fitno <<- fitno+1
     return(sim.chisq)
  }
  
  # fit glm to data
  binomial.mle <- glm.fit(object$X, object$Y,
                          offset = object$offset, family = binomial())
  
  
  thedata <- list(X=object$X,Y=object$Y,offset=object$offset,mle=binomial.mle)
  fitno <- 1
  binomial.boot <- boot(thedata, binomial.fun, R, sim="parametric",
                       ran.gen = binomial.rg, mle=binomial.mle, parallel=parallel,
                       ncpus = cores)
  return(binomial.boot)
}

outlierTest <-
  ## Short form for generic function
  function(object, R = 999, cores = max(detectCores() - 1, 1))
    UseMethod("outlierTest")

outlierTest.robmixglm <- function(object, R = 999, cores = max(detectCores() - 1, 1)) {
  if (!inherits(object, "robmixglm"))
    stop("Use only with 'robmixglm' objects.\n")
  if (cores>1) {
    if(.Platform$OS.type=="unix") parallel <- "multicore"
    else parallel <- "snow"
  } else parallel <- "no"
  out <- switch (
    object$family,
    gaussian = gaussian.outlierTest.robmixglm(object, R, parallel, cores),
    binomial = binomial.outlierTest.robmixglm(object, R, parallel, cores),
    poisson = poisson.outlierTest.robmixglm(object, R, parallel, cores),
    gamma = gamma.outlierTest.robmixglm(object, R, parallel, cores),
    truncpoisson = truncpoisson.outlierTest.robmixglm(object, R, parallel, cores),
    nbinom = nbinom.outlierTest.robmixglm(object, R, parallel, cores)
  )
  class(out) <- "outlierTest"
  return(out)
}
