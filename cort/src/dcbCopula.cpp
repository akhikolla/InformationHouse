#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::NumericVector dcbCopula(const NumericMatrix u,
                              const NumericMatrix x,
                              const NumericVector m){

  int n_obs = u.nrow();
  int dim = x.ncol();
  int n_pts = x.nrow();
  Rcpp::NumericMatrix floor_x(n_pts,dim);
  Rcpp::NumericMatrix floor_u(n_obs,dim);
  Rcpp::NumericVector rez(n_obs);
  rez.fill(0.0);

  // coerce grid :
  for(int i = 0; i < n_pts; i++){
    floor_x(i,_) = 1.0*floor(x(i,_)*m);
  }

  for (int j = 0; j < n_obs; j++){
    // coerce points to inf_points :
    floor_u(j,_) = 1.0*floor(u(j,_)*m);
    for(int i = 0; i < n_pts; i++){
      if(is_true(all(floor_u(j,_) == floor_x(i,_)))){
        rez(j) += 1.0/n_pts;
      }
    }
  }
  for(int d =0; d < dim; d++){
    rez = 1.0 * rez * m(d);
  }
  return(rez);
}
