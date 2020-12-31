#include "cpphutils.h"

// [[Rcpp::export]]
DoubleVector do_summary3_dbl(DoubleVector x, DoubleVector y, DoubleVector z,
                             bool in_place = false,
                             bool do_max = true) {
  const R_xlen_t n = x.length();
  const int ny = y.length();
  const int nz = z.length();
  const bool y1 = ny == 1;
  const bool z1 = nz == 1;
  NumericVector out = in_place ? NumericVector(x) : NumericVector(clone(x));
  const double y0 = y[0];
  const double z0 = z[0];

  if (do_max) {
    for (R_xlen_t i = 0; i < n; ++i) {
      double xi = x[i];
      double yi = y1 ? y0 : y[i];
      double zi = z1 ? z0 : z[i];
      if (xi < zi && yi < zi){
        out[i] = zi;
      } else {
        if (xi < yi){
          out[i] = yi;
        }
      }
    }
  } else {
    for (R_xlen_t i = 0; i < n; ++i) {
      double xi = x[i];
      double yi = y1 ? y0 : y[i];
      double zi = z1 ? z0 : z[i];
      if (xi > zi && yi > zi){
        out[i] = zi;
      } else {
        if (xi > yi){
          out[i] = yi;
        }
      }
    }
  }
  return out;
}


// [[Rcpp::export]]
IntegerVector do_summary3_int(IntegerVector x, IntegerVector y, IntegerVector z,
                       bool in_place = false,
                       bool do_max = true) {
  const R_xlen_t n = x.length();
  const int ny = y.length();
  const int nz = z.length();
  const bool y1 = ny == 1;
  const bool z1 = nz == 1;
  IntegerVector out = in_place ? IntegerVector(x) : IntegerVector(clone(x));
  const int y0 = y[0];
  const int z0 = z[0];

  if (do_max) {
    for (R_xlen_t i = 0; i < n; ++i) {
      int xi = x[i];
      int yi = y1 ? y0 : y[i];
      int zi = z1 ? z0 : z[i];
      if (xi < zi && yi < zi){
        out[i] = zi;
      } else {
        if (xi < yi){
          out[i] = yi;
        }
      }
    }
  } else {
    for (R_xlen_t i = 0; i < n; ++i) {
      int xi = x[i];
      int yi = y1 ? y0 : y[i];
      int zi = z1 ? z0 : z[i];
      if (xi > zi && yi > zi){
        out[i] = zi;
      } else {
        if (xi > yi){
          out[i] = yi;
        }
      }
    }
  }
  return out;
}
