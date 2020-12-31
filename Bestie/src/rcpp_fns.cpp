
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector collectC(IntegerVector xs, NumericVector ys, int n) {
  int ln = ys.size();
  NumericVector out(n);
  
  for(int i = 0; i < ln; ++i) {
    out[xs[i]] += ys[i];
  }
  return out;
}
