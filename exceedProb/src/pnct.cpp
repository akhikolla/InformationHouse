// [[Rcpp::depends(BH)]]

#include <Rcpp.h>
#include <boost/math/distributions/non_central_t.hpp> 

using namespace boost::math;

//' t-distribution with Boost
//'
//' This function returns the cdf of a noncentral t-distribution.
//' It is more accurate than stats::pt() for large ncp
//' @param x Test statistic
//' @param df Degrees of freedom
//' @param ncp Noncentrality parameter
//' @return Cumulative probability 
//' @export
// [[Rcpp::export]]
double pnct(const double x, const double df, const double ncp) {
  non_central_t dist(df, ncp);
  return cdf(dist, x);
}
