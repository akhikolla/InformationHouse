#include <Rcpp.h>
using namespace Rcpp;

//' @rdname LossMat
//' @title C++ Function for Interaction Loss Function
//' @name LossMat
//' @description This function calculates the loss function for the interaction clustering
//' for all data slices and clusters means. The inputs are numeric matrices.
//' @param x The data matrix, with the N slices strung out as vectors in the columns.
//' @param y The matrix of cluster means, with each mean represented by a row.
//' @return A numeric matrix with \code{nclust} rows and \code{N} columns.

// [[Rcpp::export]]
NumericMatrix LossMat(NumericMatrix x, NumericMatrix y) {
  // x (datmat) is J*K by N, y (mnsmat) is J*K by nclust
  int JK = x.nrow(), N = x.ncol(), nclust = y.ncol();
  NumericMatrix out(nclust, N);
  
  for (int j = 0; j < N; j++) {
  // j is over N
    for (int k = 0; k < nclust; k++) {
    // k is over nclust
      for (int i = 0; i < JK; i++) {
      // i is over JK
        out(k, j) += pow(x(i, j) - y(i, k), 2.0);
      }
    }
  }
  return out;
}
