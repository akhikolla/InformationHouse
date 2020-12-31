// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

//' Fast addition of vector to each row of matrix
//'
//' @description Fast addition of vector to each row of a matrix. This corresponds to t(t(x) + v)
//' @param x A matrix with dimensions n*k.
//' @param v A vector of length k.
//' @return A matrix of dimension n*k where v is added to each row of x
//' @author Claus Ekstrom <claus@@rprimer.dk>
//' @examples
//'
//' A <- matrix(1:12, ncol=3)
//' B <- c(1, 2, 3)
//'
//' add_torows(A, B)
//'
//' @export
// [[Rcpp::export]]
arma::mat add_torows(const arma::mat &x, const arma::rowvec &v) {

  if (x.n_cols != v.n_elem) {
     stop("Non-conformable objects. Length of v does not match columns in x");
  }

  arma::mat m = x.each_row() + v;
  
  return (m);
}


