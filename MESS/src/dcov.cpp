// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// #include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;


//' Fast distance covariance matrix
//'
//' @description Fast computation of the distance covariance between two matrices with the same number of rows.
//' @param x A matrix with dimensions n*k.
//' @param y A matrix with dimensions n*l.
//' @return A number representing the distance covariance between x and y
//' @author Claus Ekstrom <claus@@rprimer.dk>
//' @export
// [[Rcpp::export]]
double dCov(NumericMatrix x, NumericMatrix y) {

  // Check the same number of rows

  arma::mat X(x.begin(), x.nrow(), x.ncol(), false);
  arma::mat Y(y.begin(), y.nrow(), y.ncol(), false);
  
  int N = x.nrow();

  if (y.nrow() != N) {
    Rcpp::stop("the two matrices must have the same number of rows");
  }

  arma::mat G = X * X.t();
  arma::mat res(N, N, arma::fill::zeros);

  arma::mat G2 = Y * Y.t();
  arma::mat res2(N, N, arma::fill::zeros);

  for (int i = 0 ; i<N; i++) {
    for (int j = i+1 ; j<N; j++) {
      res(i, j) = sqrt(G(i,i) - 2*G(i, j) + G(j, j));
      res(j, i) = res(i,j);

      res2(i, j) = sqrt(G2(i,i) - 2*G2(i, j) + G2(j, j));
      res2(j, i) = res2(i,j);
    }
  }

  // Compute Euclidian matrix for X
  res.each_col() -= mean(res, 1);
  res.each_row() -= mean(res, 0);
  res += mean(mean(res));  // Could be optimized


  res2.each_col() -= mean(res2, 1);
  res2.each_row() -= mean(res2, 0);
  res2 += mean(mean(res2));  // Could be optimized


  return(  sqrt(accu(res % res2))/N );
}




//' Fast distance correlation matrix
//'
//' @description Fast computation of the distance correation matrix between two matrices with the same number of rows. Note that this is not the same as the correlation matrix distance that can be computed with the cmd function.
//' @param x A matrix with dimensions n*k.
//' @param y A matrix with dimensions n*l.
//' @return A number between 0 and 1 representing the distance covariance between x and y
//' @author Claus Ekstrom <claus@@rprimer.dk>
//' @export
// [[Rcpp::export]]
double dCor(NumericMatrix x, NumericMatrix y) {
  return(dCov(x, y)/sqrt(dCov(x, x)*dCov(y, y)));
}

