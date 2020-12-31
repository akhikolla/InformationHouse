#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector ransamf(int repnum, int n, NumericVector B1inq, NumericVector xi1d, NumericVector xi2d, NumericVector cids1, NumericVector cids2) {

  Rcpp::NumericVector out(repnum);


  for (int i = 0; i < repnum; ++i) {
    double sxi1dg = 0;
    double sxi2dg = 0;
    double scidsg1 = 0;
    double scidsg2 = 0;
    
    NumericVector g = rnorm(n);
  
    for(int h = 0; h < n; ++h) {
      sxi1dg += xi1d[h]*g[h];
      sxi2dg += xi2d[h]*g[h];
      scidsg1 += cids1[h]*g[h];
      scidsg2 += cids2[h]*g[h];
    }

    out[i] = (B1inq[0]*sxi1dg + B1inq[1]*sxi2dg)/std::sqrt(n) + scidsg1 + scidsg2;

  }

  return out;

}
