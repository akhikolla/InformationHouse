#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]

NumericMatrix support() {
  int nr = 4950, nc = 2, k = 0;
  double y = 0;
  NumericMatrix mat(nr,nc);
  for(int i = 0; i < 99; i++) {
     y = 0.001 + 0.01*i;
     for(int j = 0; j < 99 - i; j++) {
       mat(k, 0) = 0.001 + 0.01*j;
       mat(k, 1) = y;
       k = k + 1;
     }
  }
  return mat;
}

