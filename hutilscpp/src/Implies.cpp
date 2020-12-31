#include "cpphutils.h"

// [[Rcpp::export]]
LogicalVector Implies(LogicalVector x, LogicalVector y, bool anyNAx = true, bool anyNAy = true) {
  R_xlen_t N = x.size();
  R_xlen_t yN = y.size();
  if (N != yN) {
    stop("lengths of x and y differ.");
  }

  LogicalVector out(clone(y));
  if (!anyNAx && !anyNAy) {
    for (R_xlen_t i = 0; i < N; ++i) {
      if (!x[i]) {
        out[i] = true;
      }
    }
    return out;
  }

  LogicalVector nax(N);
  LogicalVector nay(N);
  nax = is_na(x);
  nay = is_na(y);
  for (R_xlen_t i = 0; i < N; ++i) {
    if (nax[i]) {
      if (nay[i] || !y[i]) {
        out[i] = NA_LOGICAL;
      } else {
        out[i] = true;
      }
      continue;
    }
    if (nay[i]) {
      if (!x[i]) {
        out[i] = true;
      }
      continue;
    }
    // No NAs
    if (!x[i]) {
      out[i] = true;
    }
  }
  return out;
}
