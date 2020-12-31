#include <Rcpp.h>
using namespace Rcpp;

// applies to ops 1-6 only
inline bool do_one_op_1_6(const int &op, const int &xi, const int &yi) {
  switch(op) {
  case 1:
    return xi != yi;
  case 2:
    return xi == yi;
  case 3:
    return xi >= yi;
  case 4:
    return xi <= yi;
  case 5:
    return xi >  yi;
  case 6:
    return xi <  yi;
  }
  return false; // # nocov
}

// [[Rcpp::export(rng = false)]]
IntegerVector do_whichs_16(int op, IntegerVector x, IntegerVector y, int nThread = 1) {
  if (op < 1 || op > 6) {
    stop("Internal error(do_whichs_16): op < 1 || op > 6"); // # nocov
  }
  R_xlen_t xn = x.length();
  R_xlen_t yn = y.length();

  if (xn >= INT_MAX || yn >= INT_MAX) {
    stop("Internal error(do_whichs_16): do_whichs1 only accepts short vectors."); // # nocov
  }

  const bool yn1 = yn == 1;
  const bool xye = yn == xn;
  int n = (xn >= yn) ? xn : yn;

  if (xn != yn && !yn1) {
    stop("Internal error(do_whichs_16): xn != yn && yn != 1"); // # nocov
  }

  IntegerVector out = no_init(n);
  int j = 0;
  for (int i = 0; i < n; ++i) {
    const int xi = x[i];
    const int yi = xye ? y[i] : y[0];
    bool resi = do_one_op_1_6(op, xi, yi);
    out[j] = i + 1;
    j += resi;
  }
  if (j == 0) {
    return IntegerVector(0);
  }
  IntegerVector o = no_init(j);
  for (int i = 0; i < j; ++i) {
    o[i] = out[i];
  }

  return o;
}
