#include<Rcpp.h>
#include<math.h>
using namespace Rcpp;
// [[Rcpp::export]]
double GoFHe(NumericVector x,double a,NumericVector a1) { 
  
 
  double s = 0;
  double n = x.size();
  double  e = n/sum(x);
  NumericVector y = e*x;
  std::sort(y.begin(),y.end());
  for (int j=0;j<n; j++){
    s=s-2*(exp(y[j]+a))*a1[j];
    for (int i=0;i<n; i++){s=s+1/n*1/(y[j]+y[i]+a);}
      } 
  return s;    
    
}
