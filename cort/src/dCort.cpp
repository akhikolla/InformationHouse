#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::NumericVector dCort(const NumericMatrix a,
                          const NumericMatrix b,
                          const NumericVector kern,
                          const NumericMatrix u){

  int n_pts = u.nrow();
  int dim = u.ncol();
  int n_leaves = a.nrow();
  bool is_in;
  Rcpp::NumericVector rez(n_pts);
  rez.fill(0.0);

  for (int i = 0; i < n_pts; i++){
    for (int l = 0; l < n_leaves; l++){
      is_in = true;
      for (int d = 0; d < dim; d++){
        is_in = is_in && (a(l,d) <= u(i,d)) && (u(i,d) < b(l,d));
        if(!is_in){
          break;
        }
      }
      if(is_in){
        rez(i) = kern(l);
        break;
      }
    }
  }
  return(rez);
}
