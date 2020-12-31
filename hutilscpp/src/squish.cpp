#include "cpphutils.h"

// [[Rcpp::export]]
NumericVector squishn(NumericVector x, double a, double b, bool in_place = false) {
  R_xlen_t N = x.size();
  NumericVector out = in_place ? NumericVector(x) : NumericVector(clone(x));
  for (R_xlen_t i = 0; i < N; ++i) {
    if (x[i] < a) {
      out[i] = a;
    } else if (x[i] > b) {
      out[i] = b;
    }
  }
  return out;
}

// [[Rcpp::export]]
IntegerVector squishi(IntegerVector x, int a, int b, bool in_place = false) {
  R_xlen_t N = x.size();
  IntegerVector out = in_place ? IntegerVector(x) : IntegerVector(clone(x));
  for (R_xlen_t i = 0; i < N; ++i) {
    if (x[i] < a) {
      out[i] = a;
    } else if (x[i] > b) {
      out[i] = b;
    }
  }
  return out;
}


