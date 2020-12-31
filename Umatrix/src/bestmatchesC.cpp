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
NumericVector bestmatchesC(NumericMatrix WeightVectors, NumericMatrix DataPoints, int Columns) {
  int nrow = WeightVectors.nrow();
  int ncol = WeightVectors.ncol();

  int nrowD = DataPoints.nrow();
//int ncolD = DataPoints.ncol();

  NumericMatrix result(nrowD,2); //posy(row), posx(col)
  
  // iteration over every DataPoint
  for(int k = 0; k < nrowD; k++){
    double dist=0;
    int col = 0;
    int row = 0;

    // calculate the distance to the first weightvector in the matrix
    for(int j = 0; j < ncol; j++){
      dist += pow((WeightVectors(0,j) - DataPoints(k,j)),2);
    }

    // initialize the current best match as the first weightvector in the matrix
    double smallestDistance = dist;
    int currentBestmatch=0;

    // iterate over every weightvector within the matrix
    for(int i = 1; i < nrow; i++) {
      dist = 0;
      // calculate the distance to the current weightvector
      for(int j = 0; j < ncol; j++){
        dist += pow((WeightVectors(i,j) - DataPoints(k,j)),2);
      }
      // when the current distance is smaller than the smallest distance up to this point,
      // change bestmatch to current index
      if(smallestDistance > dist){
        smallestDistance = dist;
        currentBestmatch=i;
      }
    }

    row = ((int)(currentBestmatch) / Columns) + 1;
    col = ((currentBestmatch) % Columns) + 1;
    result(k,0) = row;
    result(k,1) = col;
  }

  // add 1 to index because indices in R start with 1 instead of 0
  return(result);
}

