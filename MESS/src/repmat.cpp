// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

//' Fast replication of a matrix
//'
//' @description Fast generation of a matrix by replicating a matrix row- and column-wise in a block-like fashion
//' @param x A matrix with dimensions r*c.
//' @param nrow An integer giving the number of times the matrix is replicated row-wise
//' @param ncol An integer giving the number of times the matrix is replicated column-wise
//' @return A matrix with dimensions (r*nrow) x (c*ncol)
//' @author Claus Ekstrom <claus@@rprimer.dk>
//' @examples
//'
//' m <- matrix(1:6, ncol=3)
//' repmat(m, 2)     # Stack two copies of m on top of each other
//' repmat(m, 2, 3)  # Replicate m with two copies on top and three copies side-by-side 
//'
//' @export
// [[Rcpp::export]]
NumericMatrix repmat(NumericMatrix x, int nrow=1, int ncol=1) {

  arma::mat X(x.begin(), x.nrow(), x.ncol(), false);

  // Sanity checks
  if (nrow<1)
    Rcpp::stop("Number of rows to replicate must be >= 1");

  if (ncol<1)
    Rcpp::stop("Number of columns to replicate must be >= 1");
    
    
  return wrap(repmat( X, nrow, ncol));
}

