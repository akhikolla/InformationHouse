#include<Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericVector GoFKl(NumericVector x,double a) { 
  
  NumericVector res(2);
  double b=0;
  double c=0;
  double n = x.size();
  double  e = n/sum(x);    
  NumericVector y = e*x;
  std::sort(y.begin(),y.end());
  for (int i=0; i<(n-1);i++){ 
    b=b+sum(pow(y[i],2)*y[Range((i+1),(n-1))]);  
    c=c+sum(exp(-a*y[i])*(a*(y[Range((i+1),(n-1))]-y[i])-2));
  }
  res[0]=b;
  res[1]=c;
  return res;
}
