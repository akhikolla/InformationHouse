#define ARMA_32BIT_WORD 1
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
NumericVector ruler(NumericMatrix vR, NumericVector uR, NumericMatrix ciR) {
  int n = vR.nrow(), k = vR.ncol();
  NumericVector uM = rep(uR, n);
  arma::mat v(vR.begin(), n, k, false);
  arma::mat u(uM.begin(), k, n, false);
  arma::mat ci(ciR.begin(), k, k, false);
  arma::mat uv = arma::trans(v) - u;
  arma::colvec distances = arma::diagvec(arma::trans(uv) * ci * uv);
  return Rcpp::wrap(distances);
}

