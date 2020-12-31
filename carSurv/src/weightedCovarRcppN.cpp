#include <Rcpp.h>
using namespace Rcpp;

//' @title Estimate weighted covariance
//' @description Efficient C implementation of the sample covariance estimator. 
//' The denominator is defined as the sample size.
//' @param x Covariate without weighting (numeric vector).
//' @param y Response. The mean of the response contains weights (numeric vector).
//' @param w Weights for averaging (numeric vector).
//' @return Weighted variance (numeric scalar).
//' @author Thomas Welchowski
//' @keywords survival
//' @note There are no safety checks of input arguments. 
//' @examples 
//' # Simulate two random vectors
//' set.seed(3975)
//' x <- rnorm(100)
//' set.seed(-3975)
//' y <- rnorm(100)
//' # Calculate variance with standard R function
//' # Rescaling ensures that both calculations use same denominator "n"
//' covarEst <- cov(x, y) * (100-1) / 100
//' # Calculate weighted variance with equal weights
//' equalW <- rep(1, 100)
//' weightCovarEst <- weightedCovarRcppN(x=x, y=y, w=equalW)
//' # Output comparison
//' all.equal(covarEst, weightCovarEst)
//' # Runtime comparison
//' library(microbenchmark)
//' microbenchmark(Default=cov(x, y), New=weightedCovarRcpp(x=x, y=y, w=equalW), times=25)
//' # -> New method is multiple times faster
//' @export
// [[Rcpp::export]]
double weightedCovarRcppN(Rcpp::NumericVector x, Rcpp::NumericVector y, 
                         Rcpp::NumericVector w) {
  
  // 1. Calculate sum of weighted responses, weights and squared weights
  int unsigned xLength = x.size();
  double sumVecX = 0;
  double sumVecY = 0;
  double sumWeights = 0;
  for(int unsigned i = 0; i < xLength; i++) {
    sumVecX = sumVecX + x[i];
    sumVecY = sumVecY + y[i] * w[i];
    sumWeights = sumWeights + w[i];
  }
  
  // 2. Calculate weighted mean of X and Y
  double sampleSize = x.size();
  double meanX = sumVecX / sampleSize;
  double wMeanY = sumVecY / sampleSize;
  
  // 3. Calculate weighted variance
  double sumTemp = 0;
  for(int unsigned i = 0; i < xLength; i++) {
    sumTemp = sumTemp + w[i] * (x[i] - meanX) * (y[i] - wMeanY);
  }

  // Output
  return sumTemp / sampleSize;
}
