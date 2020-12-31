// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

//' Apply cumsum to each column of matrix
//' 
//' Fast computation of apply(m, 2, cumsum)
//'
//' @param m A matrix
//' @return A matrix the same size as m with the column-wise cumulative sums.
//' @author Claus Ekstrom <claus@@rprimer.dk>
//' @examples
//'   # Generate a 100 by 10000 matrix
//'   x <- matrix(rnorm(100*10000), nrow=100)
//'   result <- colCumSum(x)
//'
//' @export
// [[Rcpp::export]]
NumericMatrix colCumSum(NumericMatrix m) {
  arma::mat M(m.begin(), m.nrow(), m.ncol(), false);
  arma::mat result;

  result = arma::cumsum(M, 0);

  return wrap(result);
}
