#include <Rcpp.h>
using namespace Rcpp;

//[[Rcpp::export]]
NumericVector transformation(NumericVector x, NumericVector lb, NumericVector ub)
{
  int n = x.size();
  NumericVector transformed_x(n);
  
  for (int i = 0; i < n; ++i)
    transformed_x[i] = x[i] * (ub[i] - lb[i]) + lb[i];
  
  return transformed_x;
}
