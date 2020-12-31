#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::NumericVector pCort(const NumericMatrix a,
                          const NumericMatrix b,
                          const NumericVector p,
                          const NumericMatrix u){

  int n_pts = u.nrow();
  int dim = u.ncol();
  int n_leaves = a.nrow();
  double measure_in;
  Rcpp::NumericVector rez(n_pts);
  rez.fill(0.0);

  for (int i = 0; i < n_pts; i++){
    for (int l = 0; l < n_leaves; l++){
      measure_in = 1.0;
      for (int d = 0; d < dim; d++){
        measure_in *= std::min(std::max(u(i,d) - a(l,d),0.0)/(b(l,d)-a(l,d)),1.0);
        if(measure_in == 0){
          break;
        }
      }
      rez(i) += measure_in*p(l);
    }
  }
  return(rez);
}
