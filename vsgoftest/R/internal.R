#### Heading ####

## File: internal.R
## Desc: internal functions of the package maxengoftest
## Date: 2017-12-21
## R version: 3.4.3
## Required packages: fitdistrplus, Rcpp
## Script designed with and optimized for RStudio

#### Optimal estimation of differential Shannon entropy for a specified model ####

## likelihood (function): a wrapper for computing the likelihood of the sample, depending on the null distribution the test is intended to be applied to
likelihood=function(x, densfun, param){
  res <- switch(densfun,
		"dunif"=dunif(x,min=param[1],max=param[2]),
    "dnorm"=dnorm(x, mean=param[1], sd=param[2]),
		"dlnorm"=dlnorm(x,meanlog=param[1],sdlog=param[2]), 
    "dexp"=dexp(x, rate=param[1]),
    "dgamma"=dgamma(x, shape=param[1], rate=param[2]),
    "dweibull"=dweibull(x, shape=param[1], scale=param[2]),
    "dpareto"=dpareto(x, param[1], param[2]),
    "df"=df(x,df1=param[1],df2=param[2],ncp=param[3]),
    "dlaplace"=dlaplace(x,mu=param[1],b=param[2]),
		"dbeta"=dbeta(x,shape1=param[1],shape2=param[2])
  )
  return(res)
}

## vs.estimate (function): Vasicek estimate (with the optimal window) of differential Shannon entropy from a numeric sample.
vs.estimate <- function(x,densfun = NULL, param, extend, delta, relax, suppress.error = FALSE){
  n <- length(x) #sample length
  TIES=FALSE
  if ((length(unique(x)) < n)){
    warning('Ties should not be present for Vasicek-Song test')
    TIES=TRUE
  }
  xord <- sort(x) #ordered statistics
  NTmax <- max(table(x)) #Max number of ties in the sample
  if (is.null(delta)) {
    delta <- switch(densfun,
                    'dunif'=1/12,
                    'dnorm'=1/12,
                    'dlnorm'=1/12,
                    'dexp'=1/12,
                    'dgamma'=2/15,
                    'dweibull'=2/15,
                    'dpareto'=1/12,
                    'df'=2/15,
                    'dlaplace'=1/12,
                    'dbeta'=2/15)
  }
  if (extend) {bound = floor(n/2)} else {
    bound <- min(max(1,floor(n^(1/3 - delta))), n/2)
  }
  if (!TIES){
    V <- vestimates(xord,1,bound)
  }
  else{
    if(bound <= NTmax -1) 
    {Vopt <- -Inf
     mopt <- Inf
     stop('Too many ties to compute Vasicek estimate.')
    }
    else{
      V <- vestimates(xord,NTmax,bound)
    }
  }
  if (relax) {#If relax = TRUE, the constraint Vmn < empirical entropy is avoided
    res <- V
  } else {#If relax = FALSE, the constraint Vmn < empirical entropy is active
    res <- V[V< -mean(log(likelihood(x, densfun, param)))] #Keep only estimates smaller than empirical log-likelihood
  }
  if (length(res)==0 & !suppress.error) {
    stop("The sample entropy is greater than empirical maximal entropy for all possible window sizes; the sample may be too small or is unlikely to be drawn from the null distribution.")
  } else{}
  if (length(res)==0 & suppress.error) {
    Vopt <- NA
    mopt <- NA
  }
  else{
    Vopt <- max(res) #Final estimate for V
    mopt <- min(which(V==Vopt)) + NTmax - 1#optimal window size m
  }
  return(list(estimate=Vopt, window=mopt, ties=TIES))
}


#### Monte-Carlo simulation of the distribution of Vasicek estimate under normality hypothesis ####

## simulate.vasicek.dist (function): simulate a sample of Vasicek test statistic for GOF to a given parametric distribution (for later use in Monte-Carlo methods).
simulate.vs.dist <- function(n,densfun,param,B, extend, delta, relax){
  tmp <- function(){
    ech <- switch(densfun,
		              "dunif"=runif(n,min=param[1],max=param[2]),
                  "dnorm"=rnorm(n, mean=param[1], sd=param[2]),
	                "dlnorm"=rlnorm(n,meanlog=param[1],sdlog=param[2]),
                  "dexp"=rexp(n, rate=param),
                  "dgamma"=rgamma(n, shape=param[1], rate=param[2]),
                  "dweibull"=rweibull(n, shape=param[1], scale=param[2]),
                  "dpareto"=rpareto(n,param[1],param[2]),
                  "df"=rf(n,df1=param[1],df2=param[2],ncp=param[3]),
                  "dlaplace"=rlaplace(n,mu=param[1],b=param[2]),
		              "dbeta"=rbeta(n,shape1=param[1],shape2=param[2]))
    est <- - vs.estimate(ech, densfun, param, extend, delta, relax, suppress.error = TRUE)$estimate
    loglik <- mean(log(likelihood(ech, densfun, param)))
    res <- est -loglik
  }
  return(replicate(B, tmp()))
}


#### MLE of the parameters of the distribution ####

## MLE.fisher (function): returns maximum-likelihood estimates of the parameters for the Fisher distribution, with specified starting values.
MLE.fisher <- function(x){
  dens <- density(x)
  mode <- dens$x[which(dens$y==max(dens$y))]
  if((mean(x)>1) & (mode<(mean(x)/(2*mean(x)-1)))) {
    d1start <- -2/(((d2start+2)/d2start)*mode-1)
    d2start <- 2*mean(x)/(mean(x)-1)
    fitdist(x,"f",start=list(df1=d1start,df2=d2start,ncp=5))
  }
  else {fitdist(x,"f",start=list(df1=1,df2=1,ncp=1),lower=c(0,0,0))}
}   

## MLE.param (function): returns maximum-likelihood estimates of the parameters for the specified distribution.
MLE.param <- function(x, densfun){
  n <- length(x)
  res <- switch(densfun,
		"dunif"=c(min(x),max(x)),
    "dnorm"=c(mean(x), sqrt((n-1)*var(x)/n)),
		"dlnorm"=c(mean(log(x)),sqrt((n-1)*var(log(x))/n)),
                "dexp"=1/mean(x),
                "dgamma"=suppressWarnings(
                  tryCatch(fitdist(x,"gamma",lower=c(0.0001,0.0001))$estimate,
                           error=function(err){cat('Unable to compute MLE of the sample for Gamma distribution.\n')}
                  )
                ),
                "dweibull"=suppressWarnings(
                  tryCatch(fitdist(x,"weibull",lower=c(0.0001,0.0001))$estimate,
                           error=function(err){cat('Unable to compute MLE of the sample for Weibull distribution.\n')}
                  )
                ),
                "dpareto"=c(1/log(exp(mean(log(x)))/min(x)), min(x)),
                "df"=MLE.fisher(x)$estimate,
                "dlaplace"=c(mean(x),mean(abs(x-mean(x)))),
		"dbeta"=suppressWarnings(
                  tryCatch(fitdist(x,"beta",lower=c(0.0001,0.0001))$estimate,
			   error=function(err){cat('Unable to compute MLE of the sample for Beta distribution.\n')}))
  )
  return(res)
}

#### Miscellianous ###
## myprint (function): append some character strings to be printed when vs.test is executed
myprint <- function(x,y) {
  n <- length(x)
  tmp <- character(n)
  for (i in 1:n){
    tmp[i] <- paste(x[i], '=', y[i], sep='')
  }
  return(paste(tmp, collapse = ', '))
}
