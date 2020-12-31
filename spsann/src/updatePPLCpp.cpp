/*******************************************************************************
Update the distance matrix
x: coordinates of all points
dm: matrix of distances
idx: row and column of the distance matrix to be updated
*******************************************************************************/
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export(".updatePPLCpp")]]

NumericMatrix updatePPLCpp(NumericMatrix x, NumericMatrix dm, int idx) {
  int ncolx = x.ncol(), nrowx = x.nrow(), ncoldm = dm.ncol(), nrowdm = dm.nrow();
  NumericVector d(nrowx, 0.0000);
  NumericMatrix x2(1, ncolx);
  NumericMatrix res(nrowdm, ncoldm);
  
  // Get the data of the distance matrix
  // This is needed so that the object passed to 'dm' is not replaced in the 
  // global environment.
  for (int i = 0; i < nrowdm; i++) {
    for (int j = 0; j < ncoldm; j++) {
      res(i, j) = dm(i, j);
    }
  }
  
  // Get the coordinates of the new point
  idx -= 1;
  for (int i = 0; i < ncolx; i++) {
    x2(0, i) = x(idx, i);
  }
  
  // Calculate distances
  for (int i = 0; i < nrowx; i++) {
    for (int j = 0; j < ncolx; j++) {
      d[i] += (x[nrowx * j + i] - x2[j]) * (x[nrowx * j + i] - x2[j]);
    }
    d[i] = pow(d[i], 0.5);
  }
  
  // replace the values in the distance matrix
  for (int i = 0; i < nrowx; i++) {
    res(i, idx) = d[i];
    res(idx, i) = d[i];
  }
  return (res);
}
/* Testing:
rm(list = ls())
Rcpp::sourceCpp('src/updatePPLCpp.cpp')
old_x <- matrix(rnorm(8), nrow = 4)
old_dm <- as.matrix(dist(old_x))
idx <- 1
new_x <- old_x
new_x[idx, ] <- rnorm(2)
#
# the new and old data_mat are different
setequal(new_x, old_x)
new_x - old_x
new_dm2 <- as.matrix(dist(new_x))
#
# the new and the old dist_mat are different
setequal(new_dm2, old_dm)
new_dm2 - old_dm
#
new_dm1 <- .updatePPLCpp(new_x, old_dm, idx)
setequal(new_dm1, old_dm)
new_dm1 - old_dm
setequal(new_dm1, new_dm2)
new_dm1 - new_dm2
*/
