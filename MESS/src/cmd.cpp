// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

//' Correlation matrix distance
//'
//' @description Computes the correlation matrix distance between two correlation matrices
//' @param x First correlation matrix
//' @param y Second correlation matrix
//' @return Returns the correlation matrix distance, which is a value between 0 and 1. The correlation matrix distance becomes
//' zero for equal correlation matrices and unity if they differ to a maximum extent.
//' @author Claus Ekstrom \email{claus@@rprimer.dk}
//' @references Herdin, M., and Czink, N., and Ozcelik, H., and Bonek, E. (2005). \emph{Correlation matrix distance, a meaningful measure for
//' evaluation of non-stationary mimo channels}. IEEE VTC.
//' @keywords univar
//' @examples
//'
//' m1 <- matrix(rep(1, 16), 4)
//' m2 <- matrix(c(1, 0, .5, .5, 0, 1, .5, .5, .5, .5, 1, .5, .5, .5, .5, 1), 4)
//' m3 <- matrix(c(1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1), 4)
//' cmd(m1, m1)
//' cmd(m1, m2)
//' cmd(m2, m3)
//'
//' @export cmd
// [[Rcpp::export]]
double cmd(NumericMatrix x, NumericMatrix y) {
  arma::mat X(x.begin(), x.nrow(), x.ncol(), false);
  arma::mat Y(y.begin(), y.nrow(), y.ncol(), false);

  if ((x.nrow() != y.nrow()) || (x.ncol() != y.ncol()))
    Rcpp::stop("the two matrices must have the same dimensions");
  
  return(1 - trace( X * Y )/ (norm(X, "fro")*norm(Y, "fro")));  
}

