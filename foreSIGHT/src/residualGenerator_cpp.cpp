#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector residualGenerator_cpp(NumericVector RN_res, double parCor1) {
  double A = parCor1;
  double B = 1;
  int n = RN_res.size();
  NumericVector res_gen(n);
  
  res_gen[0] = B*RN_res[0];
  for(int i = 1; i < n; ++i) {
    res_gen[i] = A*res_gen[i-1]+B*RN_res[i];
  }
  
  return res_gen;
}


