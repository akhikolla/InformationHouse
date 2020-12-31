# CREATING THE LOG-POSTERIOR OBJECTIVE FUNCTION
### LAST UPDATE: NA

#' Posterior Distribution for the Kriging Process
#'
#' This function finds the posterior density of a geospatial linear regression 
#' model given a point-referenced geospatial dataset and a set of parameter values. 
#' The function is useful for finding the optimum of or for sampling from the posterior 
#' distribution.
#'
#' @param tau2 Value of the nugget, or non-spatial error variance.
#' @param phi Value of the decay term, driving the level of spatial correlation.
#' @param sigma2 Value of the partial sill, or maximum spatial error variance.
#' @param beta Coefficients from linear model.
#' @param y The outcome variable that is used in the kriging model.
#' @param X The matrix of explanatory variables used in the kriging model.
#' @param east Vector of eastings for all observations.
#' @param north Vector of northings for all observations.
#' @param semivar.exp This exponent, which must be greater than 0 and less than or 
#'   equal to 2, specifies a powered exponential correlation structure for the data. 
#'   One widely used specification is setting this to 1, which yields an exponential 
#'   correlation structure. Another common specification is setting this to 2 (the 
#'   default), which yields a Gaussian correlation structure.
#' @param p.spatial.share Prior for proportion of unexplained variance that is spatial 
#'   in nature. Must be greater than 0 and less than 1. Defaults to an even split.
#' @param p.range.share Prior for the effective range term, as a proportion of the 
#'   maximum distance in the data. Users should choose the proportion of distance 
#'   at which they think the spatial correlation will become negligible. Must be 
#'   greater than 0. Values greater than 1 are permitted, but users should recognize 
#'   that this implies that meaningful spatial correlation would persist outside 
#'   of the convex hull of data. Defaults to half the maximum distance.
#' @param p.range.tol Tolerance term for setting the effective range. At the distance 
#'   where the spatial correlation drops below this term, it is judged that the 
#'   effective range has been met. Users are typically advised to leave this at 
#'   its default value of 0.05 unless they have strong reasons to choose another 
#'   level. Must be greater than 0 and less than 1.
#' @param p.beta.var Prior for the variance on zero-meaned normal priors on the 
#'   regression coefficients. Defaults to 10.
#' @param tot.var Combined variance between the nugget and partial sill. Defaults 
#'   to the variance of y. The \code{metropolis.krige} function inserts the residual 
#'   variance from a standard linear model.
#' @param local.Sigma The user is advised to ignore this option, or leave it the 
#'   value of \code{NULL}. This option is included to reduce the number of computations 
#'   required when this function is called by \code{metropolis.krige}.
#' @param max.distance The user is advised to ignore this option, or leave it the 
#'   value of \code{NULL}. This option is included to reduce the number of computations 
#'   required when this function is called by \code{metropolis.krige}.
#' 
#' @return A single number that is the posterior density of the function, which is 
#'   stored in object of class \code{matrix}.  
#'   
#' @details  This function finds the posterior density for a kriging model. It is 
#'   designed to be an internal function but is exported in the hope of it can be 
#'   useful to some users. The function utilizes information provided about the 
#'   parameters \code{tau2}, \code{phi}, \code{sigma2}, and \code{beta}. It also 
#'   utilizes the observed data \code{y}, \code{X}, \code{east}, and \code{north}. 
#'   Given a set of parameter values as well as the observed data, the function 
#'   returns the posterior density for the specified model.
#'   
#' @references 
#'   Jeff Gill. 2020. Measuring Constituency Ideology Using Bayesian Universal Kriging. 
#'   \emph{State Politics & Policy Quarterly}. \code{doi:10.1177/1532440020930197}
#'   
#' @examples
#' # Summarize Data
#' summary(ContrivedData)
#' 
#' #Initial OLS Model
#' contrived.ols<-lm(y~x.1+x.2,data=ContrivedData);summary(contrived.ols)
#' 
#' #Define Covariate Matrix
#' covariates<-cbind(1,ContrivedData$x.1,ContrivedData$x.2)
#' 
#' # Find the posterior density for the Contrived Data if all parameters were 1:
#' s.test <- krige.posterior(tau2=1,phi=1,sigma2=1,beta=rep(1,ncol(covariates)), 
#'   y=ContrivedData$y,X=covariates,east=ContrivedData$s.1,north=ContrivedData$s.2)
#' 
#' # Print posterior density
#' s.test
#' 
#' @importFrom Rcpp evalCpp
#' @export

krige.posterior<-function(tau2,phi,sigma2,beta,y,X,east,north,semivar.exp=2,
                          p.spatial.share=0.5,p.range.share=0.5,p.range.tol=0.05,
                          p.beta.var=10,tot.var=var(y),local.Sigma=NULL,max.distance=NULL){
  if(any(is.na(y))) stop("missing values in 'y' ")
	if(TRUE %in% apply(X, 2, function(x) any(is.na(x)))) stop("missing values in 'X' ")
	if(any(is.na(east))) stop("missing values in 'east' ")
	if(any(is.na(north))) stop("missing values in 'north' ")
  if(semivar.exp<=0 | semivar.exp>2) stop("semivar.exp must be greater than 0 and less than or equal to 2.")
  if(p.spatial.share<=0 | p.spatial.share>=1) stop("p.spatial.share must be between 0 and 1.")
  if(p.range.share<=0) stop("p.range.share must be greater than 0.")
  if(p.range.tol<=0 | p.range.tol>=1) stop("p.range.tol must be between 0 and 1.")
  if(p.beta.var<=0) stop("p.beta.var must be greater than 0.")
  if(is.null(local.Sigma) | is.null(max.distance)){
	  distance<-k_distmat(cbind(east,north))
	  max.distance<-max(distance)
	  Sigma<-ifelse(distance>0, sigma2*exp(-abs(phi*distance)^semivar.exp), tau2+sigma2)
  }else Sigma<-local.Sigma
  mu<-X%*%beta
  log.det.Sigma<-sum(log(k_eigenvalue(Sigma)))
  logpost<- -(1/2)*log.det.Sigma-0.5*t(y-mu)%*%k_inv(Sigma)%*%(y-mu)-
     t(beta)%*%beta/(2*p.beta.var)-log(sigma2)*(1+2*p.spatial.share)/p.spatial.share-tot.var/sigma2+
     log(tau2)*(2*p.spatial.share-3)/(1-p.spatial.share)-tot.var/tau2-
     log(1/phi)*(1+2*p.range.share)/p.range.share-max.distance*phi/((-log(p.range.tol))^(1/semivar.exp))
  return(logpost)
}
