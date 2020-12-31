#include <Rcpp.h>
using namespace Rcpp;

// Searches for the closest Vector to DataPoint within the matrix WeightVectors using the euclidean distance
// INPUT
// WeightVectors[1:m,1:n]	    WeightVectors. n weights with m variables each
// DataPoint[1:m]		    Vector with m variables
// OUTPUT
// bestmatch			    Index of the bestmatch
// author: FL
// [[Rcpp::export]]
int bestmatchC(NumericMatrix WeightVectors, NumericVector DataPoint) {
  int nrow = WeightVectors.nrow();
  int ncol = WeightVectors.ncol();
  double dist=0;

  // calculate the distance to the first weightvector in the matrix
  for(int j = 0; j < ncol; j++){
    dist += pow((WeightVectors(0,j) - DataPoint[j]),2);
  }

  // initialize the current best match as the first weightvector in the matrix
  double smallestDistance = dist;
  int currentBestmatch=0;

  // iterate over every weightvector within the matrix
  for(int i = 1; i < nrow; i++) {
    dist = 0;
    // calculate the distance to the current weightvector
    for(int j = 0; j < ncol; j++){
      dist += pow((WeightVectors(i,j) - DataPoint[j]),2);
    }
    // when the current distance is smaller than the smallest distance up to this point,
    // change bestmatch to current index
    if(smallestDistance > dist){
      smallestDistance = dist;
      currentBestmatch=i;
    }
  }

  // add 1 to index because indices in R start with 1 instead of 0
  return(currentBestmatch+1);
}
