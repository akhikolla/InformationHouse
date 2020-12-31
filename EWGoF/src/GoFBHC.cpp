#include<Rcpp.h>
#include<math.h>
using namespace Rcpp;
// [[Rcpp::export]]
double GoFBHC(NumericVector x,double a) { 
  
 
  double s = 0;
  double n = x.size();
  double  e = n/sum(x);
  NumericVector y = e*x;
  std::sort(y.begin(),y.end());
  for (int k=0;k<n; k++){
    for (int l=0;l<n; l++){
          s = s + 2 - 3*exp(-std::min(y[k],y[l])) - 2*std::min(y[k],y[l])*(exp(-y[k])+exp(-y[l]));
          s = s +2*exp(-std::max(y[k],y[l]));
      }
      } 
      s = s * 1/n;
      return s;
}
