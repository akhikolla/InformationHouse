##################################
# Estimate car scores for survival

#' @title Estimate Correlation-Adjusted Regression Survival (CARS) Scores
#' @description Estimates CARS scores. CARS scores measure the relative importance of each 
#' variable with respect to the survival times adjusted by IPC weighting.
#' @param obsTime Observed time points of a right censored survival process
#' (numeric vector).
#' @param obsEvent Observed event indicator of right censored survival process
#' (numeric vector) 0=no event, 1=event
#' @param X Data of design variables (numeric matrix). Must be already encoded.
#' @param maxIPCweight Specifies the maximum possible weight, 
#' to ensure numerical stability.
#' @param denom Specifies the denominator of the weighted sums. Two options are available: 
#' The default value "1/n" uses the sample size as denominator. Option "sum_w" uses the sum of all IPC weights in the denominator.
#' @return Estimated CAR survival score of each variable (numeric vector).
#' @details CARS scores are defined as theta=P_X^(-1/2) P_(X, log(T)). 
#' The term P_X^(-1/2) is the inverse square root of the correlation matrix between 
#' covariates X. P_(X, log(T)) is the correlation vector between covariates and 
#' the logarithmic survival time adjusted for censoring by IPC weighting.
#' @author Thomas Welchowski
#' @keywords survival
#' @note It is recommended to use default setting "denom=1/n" because in this case 
#' CARS scores are consistent. Furthermore the simulation results of "1/n" show lower 
#' root mean squared error of CARS scores with respect to the true parameter.
#' @references
#' Welchowski, T. and Zuber, V. and Schmid, M., (2018), Correlation-Adjusted Regression Survival Scores for High-Dimensional Variable Selection, <arXiv:1802.08178>
#' 
#' Zuber, V. and Strimmer, K., (2011), High-Dimensional Regression and Variable 
#' Selection Using CAR Scores, Statistical Applications in Genetics and Molecular Biology
#' 
#' Schaefer, J. and Strimmer, K., (2005), A Shrinkage Approach to Large-Scale Covariance Matrix Estimation and Implications for Functional Genomics,
#' Statistical Applications in Genetics and Molecular Biology
#' 
#' Van der Laan, M. J. and Robins, J. M., (2003), Unified Methods for Censored Longitudinal Data and Causality, Springer Series in Statistics
#' @seealso \code{\link{carVarSelect}}
#' @examples 
#' ##########################################
#' # Simulate accelerated, failure time model
#' 
#' # Generate multivariate normal distributed covariates
#' noObs <- 100
#' noCovar <- 10
#' library(mvtnorm)
#' set.seed(190)
#' X <- rmvnorm(noObs, mean=rep(0, noCovar), sigma=diag(noCovar))
#' 
#' # Generate gamma distributed survival times
#' # Only the first 5 variables have an influence
#' eta <- 1 - 2 * X[,1] - X[,2] + X[,3] +
#' 0.5 * X[,4] + 1.5 * X[,5]
#' 
#' # Function to generate survival times
#' genSurv <- function(x) {
#' set.seed(x)
#' rgamma(1, shape=2, scale=exp(eta[x]))
#' }
#' 
#' # Generate survival times
#' survT <- sapply(1:length(eta), genSurv)
#' 
#' # Generate exponential distributed censoring times
#' censT <- rexp(noObs, rate=1)
#' 
#' # Calculate event indicator
#' eventInd <- ifelse(survT <= censT, 1, 0)
#' 
#' # Calculate observed times
#' obsTime <- survT
#' obsTime[survT > censT] <- censT [survT > censT]
#' 
#' # Estimate CAR scores
#' carScores <- carSurvScore(obsTime=obsTime, obsEvent=eventInd, X=X)
#' carScores
#' @export
carSurvScore <- function(obsTime, obsEvent, X, maxIPCweight=10, denom="1/n") {
  
  # 0. Transform obsEvent to the logarithmic scale
  obsTime <- log(obsTime)

  # 1. Estimate correlations between response and covariates 
  # using IPC weighting
  
  # 1.1. Estimate columnwise variances of X
  varX <- sapply(1:dim(X)[2], function(j) var(X[, j]))
  
  # 1.2. Estimate IPC weights
  ipcWeights <- IPCweights(x=Surv(time=obsTime, event=obsEvent), 
                           maxweight=maxIPCweight)

  if(denom=="1/n") {
    # 1.3. Estimate weighted variance of time
    weightedVar <- weightedVarRcppN(y=obsTime, w=ipcWeights)
    
    # 2. Estimate weighted covariance between y and all columns of X
    weightedCovar <- sapply(1:dim(X)[2], function(j) weightedCovarRcppN(x=X[, j], y=obsTime, w=ipcWeights))
  }
  
  if(denom=="sum_w") {
    # 1.3. Estimate weighted variance of time
    weightedVar <- weightedVarRcpp(y=obsTime, w=ipcWeights)
    
    # 2. Estimate weighted covariance between y and all columns of X
    weightedCovar <- sapply(1:dim(X)[2], function(j) weightedCovarRcpp(x=X[, j], y=obsTime, w=ipcWeights))
  }
  
  # 3. Estimate correlations between y and all columns of X
  R_XY <- weightedCovar / (sqrt(weightedVar) * sqrt(varX))
  
  # 4. Compute the CAR survival scores
  CAR <- c(crossprod.powcor.shrink(x=X, y=R_XY, alpha=-0.5, verbose=FALSE))
  
  # 5. Set values outside the range of correlation coefficient to nearest number -1 or 1
  CAR[which(CAR > 1)] <- 1
  CAR[which(CAR < -1)] <- -1
  return(CAR)
}
