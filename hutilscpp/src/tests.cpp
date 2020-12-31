#include "cpphutils.h"

// [[Rcpp::export(rng = false)]]
bool test_single_ox_x1_x2(SEXP x, int op, SEXP x1, SEXP x2) {
  bool o = false;
  switch(TYPEOF(x)) {
  case INTSXP: {
    IntegerVector x_ = x;
    IntegerVector x1_ = x1;
    IntegerVector x2_ = x2;
    return single_ox_x1_x2(x_[0], op, x1_[0], x2_[0]);
  }
  case REALSXP: {
    NumericVector xd_ = x;
    NumericVector xd1_ = x1;
    NumericVector xd2_ = x2;
    return single_ox_x1_x2(xd_[0], op, xd1_[0], xd2_[0]);
  }
  }
  return o; // # nocov
}
