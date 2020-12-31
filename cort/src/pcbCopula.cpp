#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::NumericVector pcbCopula(const NumericMatrix u,
                              const NumericMatrix x,
                              const NumericVector m){

  int n_obs = u.nrow();
  int dim = x.ncol();
  int n_pts = x.nrow();
  double measure_in;
  Rcpp::NumericVector inf_point(dim);
  Rcpp::NumericMatrix seuil_inf(n_pts,dim);
  Rcpp::NumericVector rez(n_obs);
  rez.fill(0.0);

  // For each leaf :
  for(int i = 0; i < n_pts; i++){
    inf_point = 1.0*floor(x(i,_)*m)/m;

    // now for each observation,
    for(int j = 0; j < n_obs; j++){
      // check it's measure inside the leaf :
      measure_in = 1.0;
      for(int d = 0; d < dim; d++){
        measure_in *= std::min(std::max(1.0*(0.0 +u(j,d) - inf_point(d))*m(d),0.0),1.0);
        if(measure_in == 0.0){
          break;
        }
      }
      rez(j)+=measure_in/n_pts;
    }
  }
  return(rez);
}
