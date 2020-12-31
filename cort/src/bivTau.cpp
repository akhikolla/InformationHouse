#include <Rcpp.h>
using namespace Rcpp;

double B(const double ak,
         const double bk,
         const double al,
         const double bl){

  double x = std::max(ak,al);
  double y = std::min(bk,bl);
  double z = std::max(ak,bl);
  double B = 0.0;
  if(x < y){
    B += y *( y/2 - al) - x * (x/2 -al);
  }
  if(z < bk){
    B += (bk - z) * (bl - al);
  }
  return(B);
}


// [[Rcpp::export]]
Rcpp::NumericMatrix bivTau(const NumericMatrix a,
                           const NumericMatrix b,
                           const NumericVector p) {
  int dim = a.ncol();
  int n_leaves = a.nrow();
  double rez;
  Rcpp::NumericMatrix tau(dim,dim);
  Rcpp::NumericVector measures(n_leaves);
  Rcpp::NumericVector kernel(n_leaves);


  for(int i = 0; i < (dim-1); i++){
    tau(i,i) = 1.0;
    for (int j = i+1; j < dim; j++){

      measures = (b(_,i) - a(_,i))*(b(_,j) - a(_,j));
      kernel = p/measures;
      rez = 0.0;

      for(int k = 0; k < n_leaves; k++){
        if(kernel(k) != 0){
          for(int l = 0; l < n_leaves; l++){
            if(kernel(l)!= 0){
              rez += B(a(k,i),b(k,i),a(l,i),b(l,i))*B(a(k,j),b(k,j),a(l,j),b(l,j))*kernel(k)*kernel(l);
            }
          }
        }
      }
      tau(i,j) = 4.0 * rez - 1;
      tau(j,i) = tau(i,j);
    }
  }
  tau(dim-1,dim-1) = 1.0; // The loop does not go on the last point, we need to add it up.
  return(tau);
}
