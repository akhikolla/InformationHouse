#include "cpphutils.h"

// [[Rcpp::export(rng = false)]]
R_xlen_t do_which_isnt_integerish(DoubleVector x) {
  R_xlen_t N = x.length();
  for (R_xlen_t i = 0; i < N; ++i) {
    double xi = x[i];
    if (!R_finite(xi) ||
        xi < -2147483647 ||
        xi > 2147483647) {
      return i + 1;
    }
    int xiint = (int)xi;
    if (xiint != xi) {
      return i + 1;
    }
  }

  return 0;
}
