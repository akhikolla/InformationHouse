#include <Rcpp.h>
#include <math.h>

using namespace Rcpp;

// [[Rcpp::export]]
LogicalMatrix fun_sublr(NumericVector u, NumericVector l, NumericVector r) {

  int nn = l.size();
  int kk = u.size();

  LogicalMatrix result(nn, kk);

  for (int n = 0; n < nn; ++n) {

    for (int k = 0; k < kk; ++k) {

      result(n, k) = (r(n) < 100000000000000) & (u[k] > l(n)) & (u[k] <= r(n));

    }

  }

  return result;

}

