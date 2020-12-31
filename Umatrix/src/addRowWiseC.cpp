#include <Rcpp.h>
using namespace Rcpp;


// Adds the Vector DataPoint to every row of the matrix WeightVectors
// INPUT
// WeightVectors[1:m,1:n]	    WeightVectors. n weights with m components each
// DataPoint[1:m]		    Vector with m components
// OUTPUT
// WeightVectors[1:m,1:n]		    
// author: FL
//
// [[Rcpp::export]]
NumericMatrix addRowWiseC(NumericMatrix WeightVectors, NumericVector DataPoint) {
  int nrow = WeightVectors.nrow();
  int ncol = WeightVectors.ncol();

  NumericMatrix mat(nrow,ncol);

  for(int i = 0; i < nrow; ++i) {
    for(int j = 0; j < ncol; ++j) {
      mat(i,j) = WeightVectors
  (i,j) + DataPoint[j];
    }
  }

  return(mat);
}
