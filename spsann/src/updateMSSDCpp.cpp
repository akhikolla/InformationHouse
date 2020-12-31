/*******************************************************************************
Update the matrix of distances
x1: coordinates of the prediction points
x2: coordinates of the new point
dm: matrix of distances
idx: column of the matrix of distances to be updated
*******************************************************************************/
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export(".updateMSSDCpp")]]

NumericMatrix updateMSSDCpp(NumericMatrix x1, NumericMatrix x2, 
                            NumericMatrix dm, int idx) {
  int ncolx1 = x1.ncol(), nrowx1 = x1.nrow(), ncoldm = dm.ncol(), 
      nrowdm = dm.nrow();
  NumericVector d(nrowx1, 0.0000);
  NumericMatrix res(nrowdm, ncoldm);
  
  // Get the data of the distance matrix
  // This is needed so that the object passed to 'dm' is not replaced in the 
  // global environment.
  for (int i = 0; i < nrowdm; i++) {
    for (int j = 0; j < ncoldm; j++) {
      res(i, j) = dm(i, j);
    }
  }
  
  for (int i = 0; i < nrowx1; i++) {
    for (int j = 0; j < ncolx1; j++) {
      d[i] += pow(x1[nrowx1 * j + i] - x2[j], 2);
    }
    res(i, idx - 1) = pow(d[i], 0.5);
  }
  return (res);
}
/* # Testing
rm(list = ls())
Rcpp::sourceCpp('src/updateMSSDCpp.cpp')
require(SpatialTools)
x1 <- matrix(rep(1, 6), ncol = 2)
x1[2, ] <- c(2, 2); x1
b <- x1[c(1, 3), ]; b
dm <- SpatialTools::dist2(x1, b); dm
idx <- as.integer(2); idx
x2 <- matrix(x1[idx, ], nrow = 1); x2
.updateMSSDCpp(x1, x2, dm, idx)
*/
