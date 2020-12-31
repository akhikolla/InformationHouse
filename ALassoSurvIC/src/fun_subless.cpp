#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
LogicalMatrix fun_subless(NumericVector u, NumericVector lessthan) {

  int nn = lessthan.size();
  int kk = u.size();

  LogicalMatrix result(nn, kk);

  for (int n = 0; n < nn; ++n) {

    for (int k = 0; k < kk; ++k) {

      result(n, k) = (u[k] <= lessthan(n));

    }

  }

  return result;

}
