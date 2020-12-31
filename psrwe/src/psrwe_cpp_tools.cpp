#include <Rcpp.h>
#include <Rmath.h>

using namespace Rcpp;

// [[Rcpp::plugins("cpp11")]]

// [[Rcpp::init]]
void my_package_init(DllInfo *dll) {
  // initialization code here
  R_useDynamicSymbols(dll, TRUE);
}


// Test Rcpp function
//
//
// @param test test parameter
//
// @export
// [[Rcpp::export]]
double crtTest(double test) {
  return pow(test,2);
}

// Match by nearest neighbor
//
// Match subjects in group candidate with subject in the target group
//
// @param target  vector of data from the target group
// @param candidates vector of data from the candidate group
// @param ratio matching ratio
//
// [[Rcpp::export]]
NumericVector c_match(NumericVector target, NumericVector candidate, int ratio) {

  NumericVector dist(candidate.size());
  NumericVector inx(candidate.size());

  int i, j, k, start;

  // initialize index of candidates
  for (i = 0; i < inx.size(); i++) {
    inx[i] = i;
  }

  // find neighbors
  start = 0;
  for (i = 0; i < target.size(); i++) {

    // calculate distance
    for (j = start; j < inx.size(); j++) {
      dist[inx[j]] = fabs(target[i] - candidate[inx[j]]);
    }

    // sort distance
    std::sort(inx.begin() + start, inx.end(),
              [&](const int& a, const int& b) {
                return (dist[a] < dist[b]);
              });

    // keep ratio neighbors
    start += ratio;
  }

  // return
  return(inx);
}
