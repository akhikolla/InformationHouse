#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

//' Calculate AUC and its sampling standard deviation (Internal C++ implementation)
//' 
//' @param fcst numeric vector of forecasts (NAs are not allowed)
//' @param obs vector of binary observations (obs[t] evaluates to TRUE if event happens at instance t, to FALSE otherwise) 
//' @return AUC and its sampling standard deviation
//' @seealso Auc AucDiff
//' @export
// [[Rcpp::export]]
NumericVector auc_cpp(NumericVector fcst, NumericVector obs) {

  int L = obs.size();
  arma::uvec i_ord = arma::sort_index(Rcpp::as<arma::vec>(fcst), "ascend");

  double n,m,nn,mm,i,j,jp1;
  double sumV, sumV2, sumW, sumW2, theta, v, w, sd_auc;

  sumV = sumV2 = sumW = sumW2 = 0.0;
  n = m = i = 0;

  while (1) {
    nn = mm = 0;
    while (1) {
      j = i_ord[i];
      if (obs[j]) {
        mm++;
      } else {
        nn++;
      }
      if (i == L-1) {
        break;
      } 
      jp1 = i_ord[i+1];
      if (fcst[j] != fcst[jp1]) {
        break;
      } 
      i++;
    }
    sumW += nn * (m + mm/2.0);
    sumW2 += nn * (m + mm/2.0) * (m + mm/2.0);
    sumV += mm * (n + nn/2.0);
    sumV2 += mm * (n + nn/2.0) * (n + nn/2.0);
    n += nn;
    m += mm;
    i++;
    if (i >= L) {
      break;
    }
  }

  theta = sumV / (m*n);
  v = sumV2 / ((m-1)*n*n) - sumV*sumV / (m*(m-1)*n*n);
  w = sumW2 / ((n-1)*m*m) - sumW*sumW / (n*(n-1)*m*m);

  sd_auc = sqrt(v / m + w / n);

  return NumericVector::create(theta, sd_auc);

}

