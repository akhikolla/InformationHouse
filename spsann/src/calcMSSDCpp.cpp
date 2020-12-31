/*******************************************************************************
Calculate the mean squared shortest distance
x: matrix of distances between candidate locations and sample points
*******************************************************************************/
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export(".objMSSD")]]
double objMSSDCpp(NumericMatrix x) {
  
  int nrow = x.nrow(), ncol = x.ncol();
  NumericVector out(nrow), xi(ncol);
  
  for (int i = 0; i < nrow; i++) {
    for (int j = 0; j < ncol; j++) {
      xi[j] = x(i, j);
    }
    out[i] = min(xi);
  }
  return mean(pow(out, 2));
}
