#include <Rcpp.h>
using namespace Rcpp;

//[[Rcpp::export]]
NumericVector anti_transformation(NumericVector x, NumericVector lb, NumericVector ub)
{
  int n = x.size();
  NumericVector transformed_x(n);
  
  for (int i = 0; i < n; ++i)
    transformed_x[i] = (x[i] - lb[i]) / (ub[i] - lb[i]);
  
  return transformed_x;
}
