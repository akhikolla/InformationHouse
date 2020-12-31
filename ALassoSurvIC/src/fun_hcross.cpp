#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix fun_hcross(NumericMatrix x) {

  int nn = x.nrow();
  int pp = x.ncol();
  int pp2 = pp * (pp + 1)/2;
  int ele = 0;

  NumericMatrix result(pp2, nn);

  for (int p1 = 0; p1 < pp; ++p1) {

    for (int p2 = 0; p2 <= p1; ++p2) {

      for (int n = 0; n < nn; ++n) {

        result(ele, n) = x(n, p1) * x(n, p2);

      }

      ele += 1;

    }

  }

  return result;

}
