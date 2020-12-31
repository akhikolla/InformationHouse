#include "cpphutils.h"

bool do_is_safe2int(double x) {
  return R_finite(x) && x <= 2147483647 && x >= -2147483647 && ((int)x == x);
}

int type_safe2int(double x) {
  if (ISNAN(x)) {
    return 2;
  }
  if (x < -2147483647 || x > 2147483647) {
    return 0;
  }
  int xi = (int)x;
  return (xi == x) ? 1 : 0;
}


//' @noRd
//' @param x Candidate vector.
//' @return
//'   0 if unsafe to coerce to integer
//'   1 if   safe to coerce to integer and _zero_ NAs in output
//'   2 if   safe to coerce to integer but _some_ NAs in output
// [[Rcpp::export(rng = false)]]
int is_safe2int(DoubleVector x) {
  R_xlen_t n = x.length();
  int out = 1;
  for (R_xlen_t i = 0; i < n; ++i) {
    double xi = x[i];
    // (int)NaN is UBD

    if (R_finite(xi) && xi <= 2147483647 && xi >= -2147483647) {
      int xint = (int)xi;
      if (xint != xi) {
        return 0; // integer not possible
      }
      continue;
    }
    if (R_IsNA(xi) || R_IsNaN(xi)) {
      out = 2; // out = 1 no longer possible
      continue;
    }
    if (!R_finite(xi)) {
      return 0;
    }
    if (xi > 2147483647) {
      return 0;
    } else if (xi + 2147483647 <= 0) {
      return 0;
    }
  }
  return out;
}

// [[Rcpp::export(rng = false)]]
IntegerVector force_as_integer(DoubleVector x, int na_code) {
  R_xlen_t n = x.size();
  IntegerVector out = no_init(n);
  // (int) non-finites are UBD
  switch(na_code) {
  case 0:
    stop("Internal error: na_code = 0 so cannot return IntegerVector safely.");
  case 1:
    for (R_xlen_t i = 0; i < n; ++i) {
      out[i] = (int)x[i];
    }
  case 2:
    for (R_xlen_t i = 0; i < n; ++i) {
      out[i] = (R_finite(x[i])) ? (int)x[i] : NA_INTEGER;
    }
  }
  return out;
}

