# include <Rcpp.h>
# include "fittingmeasure.h"
using namespace Rcpp;
///' Weighted sampling without replacement from a finite urn
///'
///'
///' @param n Number of distinct integer-labeled balls in the urn.
///' @param size Number of balls sampled without replacement.
///' @param prob Numeric vector of length \eqn{n} probability masses assigned to each ball.
///'
///' @return vector
// [[Rcpp::export]]
Rcpp::IntegerVector quickintsample(int n, int size, Rcpp::NumericVector prob){

  int num = size ;
  int x = n;
Rcpp::IntegerVector vx = Rcpp::clone<Rcpp::IntegerVector>(x);
Rcpp::NumericVector pr = Rcpp::clone<Rcpp::NumericVector>(prob);
Rcpp::NumericVector rnd = rexp(x) / pr;
for(int i= 0; i<vx.size(); ++i) vx[i] = i;
std::partial_sort(vx.begin(), vx.begin() + num, vx.end(), Comp(rnd));
vx = vx[seq(0, num - 1)] + 1;
return vx;
}
// Can we optimize for speed avoiding the draw of n random values in rexp(x)?

