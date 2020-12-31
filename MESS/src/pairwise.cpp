// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

//' Compute Schur products (element-wise) of all pairwise combinations of columns in matrix
//'
//' Fast computation of all pairwise element-wise column products of a matrix.
//'
//' Note that the output order of columns corresponds to the order of the columns in x. First column 1 is multiplied with each of the other columns, then column 2 with the remaining columns etc. 
//'
//' @param x A matrix with dimensions r*c.
//' @param self A logical that determines whether a column should also be multiplied by itself.
//' @return A matrix with the same number of rows as x and a number of columns corresponding to c choose 2 (+ c if self is TRUE), where c is the number of columns of x. 
//' @author Claus Ekstrom <claus@@rprimer.dk>
//' @examples
//'
//' X <- cbind(rep(1, 4), 1:4, 4:1)
//' pairwise_Schur_product(X)
//' pairwise_Schur_product(X, self=TRUE)
//'
//' @export
// [[Rcpp::export]]
NumericMatrix pairwise_Schur_product(NumericMatrix x, bool self=false) {

  arma::mat X(x.begin(), x.nrow(), x.ncol(), false);
  arma::mat res(x.nrow(), (x.ncol()*(x.ncol()-1)/2  + ((self) ? x.ncol() : 0 )));

  int index = 0;
  for (int i=0; i<(x.ncol()-1 + ((self) ? 1 : 0)); i++) {
    for (int j=(i+1 + ((self) ? -1 : 0)); j<x.ncol(); j++) {
      res.col(index) = X.col(i) % X.col(j);
      index++;
    }
  }
  
  return wrap(res);
}




//' Compute all pairwise combinations of indices
//'
//' Fast computation of indices of all pairwise element of a vector of length n.
//'
//' Note that the output order of columns corresponds to the order of the columns in x. First column 1 is multiplied with each of the other columns, then column 2 with the remaining columns etc. 
//'
//' @param n A number giving the number of elements to create all pairwise indices from
//' @param self A logical that determines whether a column should also be multiplied by itself.
//' @return A matrix with n*(n+1)/2 rows (if self=TRUE) or n*(n-1)/2 rows (if self=FALSE, the default) and two columns gicing all possible combinations of indices.
//' @author Claus Ekstrom <claus@@rprimer.dk>
//' @examples
//'
//' pairwise_combination_indices(3)
//' pairwise_combination_indices(4, self=TRUE)
//'
//' @export
// [[Rcpp::export]]
NumericMatrix pairwise_combination_indices(unsigned long n, bool self=false) {

  if (n<1)
    Rcpp::stop("The number of indices must be greater than 0");

  unsigned long nrows = (self) ? n*(n+1)/2 : n*(n-1)/2;
  
  arma::mat res(nrows, 2);

  unsigned long col1 = 1;
  unsigned long col2start = (self) ? 1 : 2;
  unsigned long col2 = col2start;
  for (unsigned long i=0; i<nrows; i++) {
    res(i, 0) = col1;
    res(i, 1) = col2;
    col2++;
    if (col2>n) {
      col2start++;
      col2=col2start;
      col1++;
    }      
  }
  
  return wrap(res);
}



/*

//' Compute Schur products (element-wise) of all pairwise combinations of columns in matrix
//'
//' Fast computation of all pairwise element-wise column products of a matrix.
//'
//' Note that the output order of columns corresponds to the order of the columns in x. First column 1 is multiplied with each of the other columns, then column 2 with the remaining columns etc. 
//'
//' @param x A matrix with dimensions r*c.
//' @param self A logical that determines whether a column should also be multiplied by itself.
//' @return A matrix with the same number of rows as x and a number of columns corresponding to c choose 2 (+ c if self is TRUE), where c is the number of columns of x. 
//' @author Claus Ekstrom <claus@@rprimer.dk>
//' @examples
//'
//' X <- cbind(rep(1, 4), 1:4, 4:1)
//' pairwise_Schur_product(X)
//' pairwise_Schur_product(X, self=TRUE)
//'
//' @export
// [[Rcpp::export]]
NumericMatrix pairwise_dcor(NumericMatrix x, bool self=false) {

  arma::mat X(x.begin(), x.nrow(), x.ncol(), false);
  arma::mat res(x.nrow(), (x.ncol()*(x.ncol()-1)/2  + ((self) ? x.ncol() : 0 )));

  /*
  arma::mat res(x.nrow(), x.nrow());

  /  *
  // Start by computing the 
  int index = 0;
  for (int i=0; i<(x.ncol()-1 + ((self) ? 1 : 0)); i++) {
    arma::mat M1();
    for (int j=(i+1 + ((self) ? -1 : 0)); j<x.ncol(); j++) {
      res.col(index) = X.col(i) % X.col(j);
      index++;
    }
  }
  
  * /
  return wrap(res);
}

*/


