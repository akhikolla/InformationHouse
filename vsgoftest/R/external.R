#### Heading ####

## File: external.R
## Desc: external (user-available) functions of the package vsgoftest
## Date: 2017-12-21
## R version: 3.4.3
## Required packages: fitdistrplus, Rcpp
## source internal.R before sourcing external.R
## Script designed with and optimized for RStudio

#### Pareto distributions: related functions ####

## dpareto (function): density function of the Pareto distribution.
dpareto <- function(x, mu, c, log=FALSE){
  ## BEGIN CHECKING ARGUMENTS VALIDITY
  if(any(mu<=0))stop("mu must be positive")
  if(any(c<=0))stop("c must be positive")
  ## END CHECKING
  tmp <- (mu*c^mu*(x)^(-1-mu))*(x>=c)
  if (log) {tmp<-log(tmp)} else{}
  return(tmp)
}


## ppareto (function): distribution function of the Pareto distribution.
ppareto <- function(q,mu,c, lower.tail=TRUE, log.p=FALSE){
  ## BEGIN CHECKING ARGUMENTS VALIDITY
  if(any(mu<=0))stop("mu must be positive")
  if(any(c<=0))stop("c must be positive") 
  ## END CHECKING
  if(lower.tail) {res <- (1-(c/q)^mu)*(q>=c)}
  else{res <- 1-(1-(c/q)^mu)(q>=c)}
  if (log.p) {res <- log(res)} else{}
  return(res)
}

##qpareto (function):  quantile function of the Pareto distribution. 
qpareto <- function(p,mu,c, lower.tail=TRUE, log.p=FALSE){
  ## BEGIN CHECKING ARGUMENTS VALIDITY
  if(any(mu<=0))stop("mu must be positive")
  if(any(c<=0))stop("c must be positive") 
  ## END CHECKING
  if(lower.tail) {res <- c*(1-p)^(-1/mu)} else {res <- c*p^(-1/mu)}
  if (log.p) {res <- log(res)} else{}
  return(res)
}

## rpareto (function): pseudo-random generator for Pareto distribution.
rpareto=function(n,mu,c){
  ## BEGIN CHECKING ARGUMENTS VALIDITY
  if(any(mu<=0))stop("mu must be positive")
  if(any(c<=0))stop("c must be positive") 
  ## END CHECKING
  u <- runif(n)
  res <- c*(1/((1-u)^(1/(mu))))
  return(res)
}

#### Laplace distributions, related functions ####

dlaplace <- function(x, mu, b, log = FALSE){
  ## BEGIN CHECKING ARGUMENTS VALIDITY
  if(!is.numeric(mu)) {stop("mu must be positive")} else{}
  if(any(b<=0)) {stop("c must be positive")} else{}
  ## END CHECKING
  tmp <- exp( - abs(x-mu)/b) /(2*b)
  if (log) {tmp<-log(tmp)} else{}
  return(tmp)
}

plaplace <- function(q,mu,b, lower.tail = TRUE, log.p = FALSE){
  ## BEGIN CHECKING ARGUMENTS VALIDITY
  if(!is.numeric(mu)) {stop("mu must be positive")} else{}
  if(any(b<=0)) {stop("c must be positive")} else{}
  ## END CHECKING
  cd <- (q < mu)
  res <- (!cd) + (-1)^(!cd)*exp((-1)^(!cd)*(q-mu)/b)/2
  if(!lower.tail) {res <- 1-res} else{}
  if (log.p) {res <- log(res)} else{}
  return(res)
}

qlaplace <- function(p,mu,b, lower.tail = TRUE, log.p = FALSE){
  ## BEGIN CHECKING ARGUMENTS VALIDITY
  if(!is.numeric(mu)) {stop("mu must be positive")} else{}
  if(any(b<=0)) {stop("c must be positive")} else{}
  ## END CHECKING
  cd <- (p > 1/2)
  res <- mu + (-1)^(cd)*b*log(2*(cd + (-1)^(cd)*p))
  if(!lower.tail) {res <- mu - res} else {}
  if (log.p) {res <- log(res)} else{}
  return(res)
}

rlaplace=function(n,mu,b){
  ## BEGIN CHECKING ARGUMENTS VALIDITY
  if(!is.numeric(mu)) {stop("mu must be positive")} else{}
  if(any(b<=0)) {stop("c must be positive")} else{}
  ## END CHECKING
  u <- runif(n)
  res <- qlaplace(p = u, mu, b)
  return(res)
}



#### Vasicek estimator of Shannon differential entropy ####

## entropy.estimate (function) : computes Vasicek estimate of differential Shannon entropy from a numeric sample
## Arguments:
##  x : a numeric vector
##  m : a numeric value specifying the window size; see details.
## Value: vasicek.estimate returns a numeric value being Vasicek estimate, with window size m of the differential entropy of the sample x.

entropy.estimate <- function(x,window){
  if (missing(x)) stop('Argument \"x\" is missing, with  no default value.') else {}
  if (!is.vector(x) | !is.numeric(x)) stop('Argument \"x\" must be a numeric vector') else{}
  if (missing(window)) stop('Argument \"window\" is missing, with no default value.') else {}
  if (!is.numeric(window) | !(length(window)==1)) stop('window must be postive and smaller than half of the sample length.') else{}
  window <- ceiling(window) # If a numeric value is given, it is tranformed into an integer.
  if (window <0 | window>length(x)/2) stop('the window must be postive and smaller than half of the sample length') else{}
  n <- length(x) #sample length
  TIES=FALSE
  if ((length(unique(x)) < n)){
    warning('Ties hinder Vasicek-Song test')
    TIES=TRUE
  }
    xord <- sort(x) #ordered statistics
  m <- window
  #Vasicek's estimator
  A <- numeric(n)
  if (!TIES){
    A[1:m] <- log(n*(xord[1:m+m]-xord[1])/(2*m))
    A[(m+1):(n-m-1)] <- log(n*(xord[(m+1):(n-m-1)+m]-xord[(m+1):(n-m-1)-m])/(2*m))
    A[(n-m):n] <- log(n*(xord[n]-xord[(n-m):n-m])/(2*m))
    V <- mean(A)
  }
  else{
    NTmax=max(table(x))
    if(m < NTmax) #if bound is too small and there are too many ties, impossible to compute Vasicek estimate.
    {stop('Too many ties for computing Vasicek estimate.')
    }
    else{
      A[1:m] <- log(n*(xord[1:m+m]-xord[1])/(2*m))
      A[(m+1):(n-m-1)] <- log(n*(xord[(m+1):(n-m-1)+m]-xord[(m+1):(n-m-1)-m])/(2*m))
      A[(n-m):n] <- log(n*(xord[n]-xord[(n-m):n-m])/(2*m))
      V <- mean(A)
    }
  }
  return(V)
}


#### Vasicek-Song GOF tests ####

## vs.test (function) : performs Vasicek-Song test
vs.test <- function(x, densfun, param = NULL, simulate.p.value = NULL, B = 5000, delta = NULL, extend = FALSE, relax = FALSE){
  ### CHECKING FOR ARGUMENT VALIDITY
  #densfun
  if(missing(densfun) | !is.character(densfun))
    stop("'densfun' must be supplied as a function or name")
  if (!(densfun %in% c("dunif","dnorm", "dlnorm", "dexp", "dgamma", "dweibull", "dpareto","df","dlaplace","dbeta"))) 
  {stop("Unsupported distribution. Please, report to help file of function \"vasicek.test\" to get the list of available distributions.")}
  else{}
  #x
  if (missing(x)) stop('Argument \"x\" is missing, with  no default value.') else{}
  if (!is.vector(x) | !is.numeric(x)) stop('Argument \"x\" must be a numeric vector') else{}
  n <- length(x) # Sample size
  Cond <- switch(densfun,
                 "dunif"=FALSE,
		 "dlnorm"=(any(x<0)),
		 "dnorm"=FALSE,
                 "dexp"=(any(x<0)),
                 "dgamma"=(any(x<0)),
                 "dweibull"=(any(x<0)),
                 "dpareto"=FALSE,
                 "df"=(any(x<0)),
                 "dlaplace"=FALSE,
		 "dbeta"=((any(x<=0))|(any(x>=1)))  
  )
  if (Cond) stop("\"x\": invalid values (not compatible with the specified distribution)")
  #param
  if (!is.null(param) & !is.numeric(param))
  {stop(paste("\"param\" must be a numeric vector specifying parameters for", densfun, sep=' '))}
  if (is.null(param)) {
    ESTIM <- suppressWarnings(
    tryCatch(MLE.param(x, densfun),
             error=function(err){stop('Unable to compute MLE for this model; see help file of fitdistrplus::fitdist for details.\n')}
    ))
    METHOD <- 'family'
  }
  else{ESTIM <- param}
  if(is.numeric(ESTIM)){
    Cond <- switch(densfun,
                   "dunif"=(length(ESTIM)!=2)|(ESTIM[2]<=ESTIM[1])|(any(x<ESTIM[1]))|(any(x>ESTIM[2])),
		               "dnorm"=(length(ESTIM)!=2)|(ESTIM[2]<=0),
		               "dlnorm"=(length(ESTIM)!=2)|(ESTIM[2]<=0),
                   "dexp"=(length(ESTIM)!=1)|(ESTIM<=0),
                   "dgamma"=(length(ESTIM)!=2)|(any(ESTIM<=0)),
                   "dweibull"=(length(ESTIM)!=2)|(any(ESTIM<=0)),
                   "dpareto"=(length(ESTIM)!=2)|(any(ESTIM<=0)),
                   "df"=(length(ESTIM)!=3)|(any(ESTIM<=0)),
                   "dlaplace"=(length(ESTIM)!=2)|(ESTIM[2]<=0),
		               "dbeta"=(length(ESTIM)!=2)|(ESTIM[1]<=0)|(ESTIM[2]<=0)
    )
    if (Cond) stop("\"param\": invalid parameter (not consistent with the specified distribution)")
    else {}
  }
  names(ESTIM) <- switch(densfun,
                         'dunif'=c('Min','Max'),
			                   'dnorm'=c('Mean','St. dev.'),
			                   'dlnorm'=c('Location','Scale'),
                         'dexp'='Rate',
                         'dgamma'=c('Shape', 'Rate'),
                         'dweibull'=c('Shape', 'Scale'),
                         'dpareto'=c('mu', 'c'),
                         'df'=c('df1','df2','ncp'),
                         'dlaplace'=c('Shape','Scale'),
			                   'dbeta'=c('Shape1','Shape2')
                         )
  #simulate.p.value
  if (!is.null(simulate.p.value) & !is.logical(simulate.p.value)){
    warning('Invalid for simulate.p.value. Reset to default: NULL.')
    simulate.p.value=NULL
  }
  if (is.null(simulate.p.value)){
    if (n < 80) {simulate.p.value=TRUE} else {simulate.p.value=FALSE}
    }
  #B
  if (!is.numeric(B) | !(length(B)==1) | B<=0) {stop("B must be a positive integer")} 
  else {
    #warning('B converted to a positive integer')
    B <- ceiling(B)
  } #Converts B into an integer
  #delta
  if (!is.null(delta) & (!(is.numeric(delta)) | !(length(delta)==1)) ) {
    warning('Argument delta must be a numeric value. Reset to default value.')
    delta <- NULL
  }
  if (is.numeric(delta)) {
    if (delta >= 1/3){
      warning('Argument delta must be a numeric value. Reset to default value.') #Bug
      delta <- NULL
    }
  }
  #extend
  if (!is.logical(extend) | ! length(extend) == 1){
    warning('Argument extend must be a logical value (TRUE or FALSE). Reset to default value.')
    extend <- FALSE
  }
  #relax
  if (!is.logical(relax) | ! length(relax) == 1){
    warning('Argument relax must be a logical value (TRUE or FALSE). Reset to default value.')
    relax <- FALSE
  }
  ### END CHECKING 
  DNAME <- deparse(substitute(x))
  NVALUE <- switch(densfun,
                   'dunif'='the uniform distribution',
		               'dnorm'='the normal distribution',
		               'dlnorm'='the log-normal distribution',
                   'dexp'='the exponential distribution',
                   'dgamma'='the gamma distribution',
                   'dweibull'='the Weibull distribution',
                   'dpareto'='the Pareto distribution',
                   'df'='the Fisher distribution',
                   'dlaplace'='the Laplace distribution',
		               'dbeta'='the Beta distribution')
  names(NVALUE) <- 'Distribution under null hypothesis'
  resVE <- vs.estimate(x, densfun, ESTIM, extend, delta, relax)
  STAT <- -resVE$estimate - mean(log(likelihood(x, densfun, ESTIM)))
  names(STAT) <- "Test statistic"
  PARAM <- resVE$window
  names(PARAM) <- 'Optimal window'
  #p-value computed by Monte Carlo simulation or asymptotic normality
    if (simulate.p.value){
    dist <- suppressWarnings(simulate.vs.dist(n, densfun, ESTIM, B, extend, delta, relax))
    if (any(is.na(dist))) {
      warning(paste('For', sum(is.na(dist)), 'simulations (over', B, '), entropy estimate is greater than empirical maximum entropy for all window sizes.'))
    }
    PVAL <- mean(dist>STAT, na.rm = TRUE)
    if (is.nan(PVAL)) { #It happens when all components of dist are NA.
      stop('Unable to compute Monte-Carlo estimate of the p-value.')
    }
  }
  else {
    bias <- log(2*PARAM) - log(n) + digamma(n+1) -digamma(2*PARAM) + 2*PARAM/n * sum(1/(1:(2*PARAM-1))) - 2*PARAM*sum(1/(1:(PARAM-1)))/n - 2*sum(((PARAM-1):1)/ (PARAM:(2*PARAM-2)))/n
    PVAL <- 1 -pnorm(  (6*PARAM*n)^(1/2)*(STAT-bias )) 
  }
  names(PVAL) <- 'p-value'
  if (is.null(param)) {
    METHOD <- paste("Vasicek-Song GOF test for",NVALUE, sep=' ')
    
    structure(list(statistic=STAT, null.value=NVALUE, parameter=PARAM, estimate=ESTIM,
                   p.value=PVAL, method=METHOD, data.name=DNAME, observed=x), class="htest")
  }
  else{
    METHOD <- paste("Vasicek-Song GOF test for",NVALUE, paste('with ', myprint(names(ESTIM), param), sep=''), 
                    sep=' ')
    structure(list(statistic=STAT, null.value=NVALUE, parameter=PARAM, 
                   p.value=PVAL, method=METHOD, data.name=DNAME, observed=x), class="htest")
    
  }
}

