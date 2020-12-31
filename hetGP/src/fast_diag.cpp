#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// Compute diag(A * B)
// [[Rcpp::export]]
NumericVector fast_diag(NumericMatrix A, NumericMatrix B){
  int nr = A.nrow();
  int nc = A.ncol();
  NumericVector di(nr);
  double mprod;

  for(int i = 0; i < nr; i++){
    mprod = 0;
    for(int j = 0; j < nc; j++){
      mprod += A(i, j) * B(j, i);
    }
    di(i) = mprod;
  }

  return(di);
}

// Compute A + diag(B) 
// Warning: modifies A!
// [[Rcpp::export]]
NumericMatrix add_diag(NumericMatrix A, NumericVector B){
  for(int i = 0; i < A.nrow(); i++){
    A(i,i) += B(i); 
  }
  return(A);
}
