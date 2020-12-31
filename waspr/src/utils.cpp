// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;
using namespace arma;


// Estimate the mode by finding the highest posterior density interval
//
// @param x a  sample from which to estimate the interval
// @param cip bandwidth for the algorithm, ranging from 0 to 1
//
// @section Details:
//   The bandwidth \code{cip} is set to 0.95 to find the 95% HPD interval and
//   to 0.1 to compute the mode.
//
// @return a scalar containing the estimate of the mode
//
// [[Rcpp::export]]

double hmode(NumericVector x, double cip) {

  int n, cil, chiv;
  double ln, M;

  n = x.size();
  NumericVector sx = clone(x);
  std::sort(sx.begin(), sx.end());

  // The number of values within the
  // (cip*100)% Confidence Interval
  cil = trunc(cip*n);

  // Will be the minimal value of the smallest interval.
  chiv = 0;

  // Size of the currently smallest interval.
  ln = sx[cil]-sx[0];

  for (int i=0; i < (n-cil); i++) {

    // If the smallest interval so far is larger than the
    // current, set the current as the new smallest interval.
    if (ln > (sx[i+cil]-sx[i])) {
      ln = (sx[i+cil]-sx[i]);
      chiv = i;
    }
  }

  M = (sx[chiv+cil]+sx[chiv])/2;

  return M;
}


// Find the highest density interval.
//
// @param x a  sample from which to estimate the interval
// @param cip bandwidth for the algorithm, ranging from 0 to 1
//
// @return a vector of length 2 containing the lower and upper bound of the HPD interval.
//
// The bandwidth \code{cip} is set to 0.95 to find the 95% HPD interval and to 0.1 to compute the mode.
//
// [[Rcpp::export]]

NumericVector hmodeci(NumericVector x, double cip) {

  int n, cil, chiv;
  double ln;

  n = x.size();
  NumericVector sx = clone(x);
  std::sort(sx.begin(), sx.end());

  // The number of values within the
  // (cip*100)% Confidence Interval
  cil = trunc(cip*n);

  // Will be the minimal value of the smallest interval.
  chiv = 0;

  // Length of the currently smallest interval.
  ln = sx[cil]-sx[0];

  for (int i=0; i < (n-cil); i++) {

    // If the smallest interval so far is larger than the
    // current, set the current as the new smallest interval.
    if (ln > (sx[i+cil]-sx[i])) {
      ln = (sx[i+cil]-sx[i]);
      chiv = i;
    }
  }

  NumericVector M(2);
  M[0] = sx[chiv];
  M[1] = sx[chiv+cil];

  return M;
}
