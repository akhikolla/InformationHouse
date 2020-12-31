// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

//' Fast quadratic form computation
//'
//' @description Fast computation of a quadratic form  \eqn{t(x) * M * x}.
//' @param x A matrix with dimensions n*k.
//' @param M A matrix with dimenions n*n. If it is to be inverted then the matrix should be symmetric and positive difinite (no check is done for this)
//' @param invertM A logical. If set to TRUE then M will be inverted before computations (defaults to FALSE)
//' @param transposex A logical. Should the matrix be transposed before computations (defaults to FALSE).
//' @return A matrix with dimensions k * k giving the quadratic form
//' @author Claus Ekstrom <claus@@rprimer.dk>
//' @export
// [[Rcpp::export]]
NumericMatrix quadform(NumericMatrix x, NumericMatrix M, bool invertM = false, bool transposex = false) {

  arma::mat X(x.begin(), x.nrow(), x.ncol(), false);
  arma::mat res;
  arma::mat A(M.begin(), M.nrow(), M.ncol(), false);
  
  if (invertM) {
    if (transposex) {
      res = X *  (arma::inv_sympd(A)*X.t());
    }
    else {
      res = X.t() *  (arma::inv_sympd(A)*X);
    }
  }
  else {
    if (transposex) {
      res = X * (A*X.t());
    }
    else {
      res = X.t() * (A*X);
    }
  }
  
  /* else {
     res = X.t() * (X);
     }
  */
  return wrap(res);
}

