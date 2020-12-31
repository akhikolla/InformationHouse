#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::NumericMatrix bivRho(const NumericMatrix a,
                           const NumericMatrix b,
                           const NumericVector p) {

  int dim = a.ncol();
  Rcpp::NumericMatrix rho(dim,dim);

  for(int i = 0; i < (dim-1); i++){
    rho(i,i) = 1.0;
    for (int j = i+1; j < dim; j++){
      rho(i,j) = 3.0*sum((b(_,i)+a(_,i)-2.0)*(b(_,j)+a(_,j)-2.0)*p) - 3.0;
      rho(j,i) = rho(i,j);
    }
  }
  rho(dim-1,dim-1) = 1.0; // The loop does not go on the last one, we ned to add it up.
  return(rho);
}

