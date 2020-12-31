#include "cpphutils.h"

// [[Rcpp::export]]
int AnyCharMatch(CharacterVector x, CharacterVector a, bool opposite = false) {
  const R_xlen_t n = x.size();
  for (R_xlen_t i = 0; i < n; ++i) {
    R_xlen_t j = (a.length() == x.length()) ? i : 0;
    if (opposite) {
      if (x[i] != a[j]) {
        return ++i;
      }
    } else {
      if (x[i] == a[j]) {
        return ++i;
      }
    }
  }
  return 0;
}

