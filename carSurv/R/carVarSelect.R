#############################################
# Variable selection with CAR survival scores

#' @title Variable selection with Correlation-Adjusted Regression Survival (CARS) Scores
#' @description Computes the CARS scores and selects significant variables. 
#' If the false non discovery rate (fndr) approach is used, significant and null variables
#' are distinguished by an a priori defined q-value. 
#' @param carSurvScores Estimated CAR survival scores of each variable (numeric vector).
#' See function \code{\link{carSurvScore}}.
#' @param method Gives the variable selection procedure. Default is "fndr", which is based
#' on the false non-discovery-rate. The other option is "threshold", which selects only 
#' variables above a given empirical quantile. 
#' @param plotDiag Should diagnostic plots of the null distribution be plotted? 
#' Default is FALSE (logical scalar).
#' @param threshold If method="threshold", then this specifies the quantile threshold of
#' the CAR survival scores. Every score above this threshold is then a significant variable. 
#' @return Index giving the significant variables of the original data (integer vector).
#' @note The quality of estimated, significant variables depends on the sample size 
#' and on the number of variables.
#' @author Thomas Welchowski
#' @keywords survival
#' @references
#' Strimmer, K., (2008), A unified approach to false discovery rate estimation, BMC Bioinformatics
#' @seealso \code{\link{carSurvScore}}
#' @examples 
#' ##########################################
#' # Simulate accelerated, failure time model
#' 
#' # Generate multivariate normal distributed covariates
#' noObs <- 100
#' noCovar <- 250
#' library(mvtnorm)
#' set.seed(7903)
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
#' # Conduct variable selection using fndr
#' carScores <- carSurvScore(obsTime=obsTime, obsEvent=eventInd, X=X)
#' selectedVar <- carVarSelect(carSurvScores=carScores)
#' selectedVar
#' 
#' # Check true positive and true negative rate
#' TPR <- mean(c(1:5) %in% selectedVar)
#' TNR <- mean(c(6:250) %in% setdiff(6:250, selectedVar))
#' perf <- TPR + TNR -1
#' perf
#' @export
carVarSelect <- function(carSurvScores, method="fndr",
                         plotDiag=FALSE, threshold=0.05) {

  if(all(is.finite(carSurvScores))) {
  
  if(method=="fndr") {
    # 1. Estimate threshold using fndr approach
    fdrScores <- fdrtool(x=carSurvScores, statistic="correlation", 
                         plot=plotDiag, verbose=FALSE)
    
    # 2. Select significant variables
    indSignif <- which(fdrScores$qval <= threshold)
  }

  if(method=="threshold") {
    # 3. Estimate cut-off value by given quantile
    estThres <- unname(quantile(carSurvScores, probs=1-threshold/2))
    
    # 4. Select variables above the threshold
    indSignif <- which(abs(carSurvScores) >= estThres)
  }
    
  # 4. Output
  return(indSignif)
  }
  else{
    return(NA)
  }
}
