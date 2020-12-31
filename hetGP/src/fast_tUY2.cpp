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

// Compute t(U) %*% Y
// [[Rcpp::export]]
NumericVector fast_tUY2(IntegerVector mult, NumericVector Y2){
  NumericVector res(mult.length());

  int idx = 0;
  int idxtmp = 0;

  for(int i = 0; i < Y2.length(); i++){
    res(idx) += Y2(i);
    idxtmp++;
    if(idxtmp == mult(idx)){
      idx++;
      idxtmp = 0;
    }

  }

  return(res);
}


