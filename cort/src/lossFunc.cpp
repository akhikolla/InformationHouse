#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double lossFunc(const NumericVector bp,
                const NumericMatrix bin_repr,
                const NumericMatrix z){

  int D = bin_repr.ncol();
  int d = bin_repr.nrow();
  int n = z.ncol();
  double n_pts, vol, loss=0;
  Rcpp::NumericVector min(d);
  Rcpp::NumericVector max(d);

  for (int n_box=0; n_box < D; n_box++){
    vol = 1.0;
    n_pts= 0.0;
    // Compute sides and volumes :
    for (int dim=0; dim < d; dim++){
      min(dim) = bp(dim)*bin_repr(dim,n_box);
      max(dim) = std::pow(bp(dim),1.0-bin_repr(dim,n_box));
      vol = vol*(max(dim) - min(dim));
    }
    if(vol > 0){
      // check if observations are inside the box :
      for(int n_obs=0; n_obs<n; n_obs++){
        if(is_true(all(min<=z(_,n_obs))) && is_true(all(z(_,n_obs)<max))){
          n_pts += 1.0;
        }
      }
      loss -=std::pow(n_pts/n,2)/vol;
    }
  }
  return(loss);
}
