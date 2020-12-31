#include <Rcpp.h>
using namespace Rcpp;

//' @title Estimate weighted variance
//' @description Efficient C implementation of the sample variance estimator. 
//' The denominator is defined as sum of all weights.
//' @param y Response. The mean of the response contains weights (numeric vector).
//' @param w Weights for averaging (numeric vector).
//' @return Weighted variance (numeric scalar).
//' @author Thomas Welchowski
//' @keywords survival
//' @note If all weights are set to 1, the denominator is identical to n.
//' There are no safety checks of input arguments.
//' @examples 
//' # Simulate a random vector
//' set.seed(3975)
//' x <- rnorm(100)
//' # Calculate variance with standard implementation
//' # Rescaling ensures that both calculations use same denominator "n"
//' varEst <- var(x) * (100-1) / 100
//' # Calculate weighted variance with equal weights
//' equalW <- rep(1/100, 100)
//' weightVarEst <- weightedVarRcpp(y=x, w=equalW)
//' # Output comparison
//' all.equal(varEst, weightVarEst)
//' # Runtime comparison
//' library(microbenchmark)
//' microbenchmark(Default=var(x), New=weightedVarRcppN(y=x, w=equalW), times=25)
//' # -> New method is multiple times faster
//' @export
// [[Rcpp::export]]
double weightedVarRcpp(Rcpp::NumericVector y, Rcpp::NumericVector w) {

  // 1. Calculate sum of weighted responses, weights and squared weights
  int unsigned yLength = y.size();
  double sumVec = 0;
  double sumWeights = 0;
  for(int unsigned i = 0; i < yLength; i++) {
   sumVec = sumVec + y[i] * w[i];
   sumWeights = sumWeights + w[i];
  }
  
  // 2. Calculate weighted mean
  double wMean = sumVec / sumWeights;
  
  // 3. Calculate weighted variance
  double sumTemp = 0;
  for(int unsigned i = 0; i < yLength; i++) {
    sumTemp = sumTemp + w[i] * (y[i] - wMean) * (y[i] - wMean);
  }

  // Output
  return sumTemp / sumWeights;
}
