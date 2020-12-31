#include "cpphutils.h"

// [[Rcpp::export]]
LogicalVector do_xor2(LogicalVector x, LogicalVector y, bool anyNAx = true, bool anyNAy = true) {
  R_xlen_t n = x.length();
  LogicalVector out = no_init(n);
  if (anyNAx && anyNAy) {
    for (R_xlen_t i = 0; i < n; ++i) {
      if (x[i] == NA_LOGICAL || y[i] == NA_LOGICAL) {
        out[i] = NA_LOGICAL;
      } else {
        bool xi = x[i];
        bool yi = y[i];
        out[i] = xi xor yi;
      }
    }
  }
  if (anyNAx && !anyNAy) {
    for (R_xlen_t i = 0; i < n; ++i) {
      if (x[i] == NA_LOGICAL) {
        out[i] = NA_LOGICAL;
      } else {
        bool xi = x[i];
        bool yi = y[i];
        out[i] = xi xor yi;
      }
    }
  }
  if (!anyNAx && anyNAy) {
    for (R_xlen_t i = 0; i < n; ++i) {
      if (y[i] == NA_LOGICAL) {
        out[i] = NA_LOGICAL;
      } else {
        bool xi = x[i];
        bool yi = y[i];
        out[i] = xi xor yi;
      }
    }

  }
  if (!anyNAx && !anyNAy) {
    for (R_xlen_t i = 0; i < n; ++i) {
      bool xi = x[i];
      bool yi = y[i];
      out[i] = xi xor yi;
    }
  }
  return out;
}




