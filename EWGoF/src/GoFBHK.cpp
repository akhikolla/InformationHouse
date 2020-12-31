#include<Rcpp.h>
#include<math.h>
using namespace Rcpp;
// [[Rcpp::export]]
double GoFBHK(NumericVector x,double a) { 
  
  double s = 0;
  double n = x.size();
  double  e = n/sum(x);    
  NumericMatrix k(2,n);
  double res=0;
  NumericVector k1(n);
  NumericVector k2(n);
  NumericVector y = e*x;
  std::sort(y.begin(),y.end());
  k(0,0) = y[1];
  k(1,0) = 0;
  for(int i=0; i<(n-1);i++){
     s = 0;
     s = sum(y[Range(0,i)]);
     k(0,(i+1)) = s/n + y[i+1]*(1-(i+1)/n) - (i+1)/n;
     k(1,(i+1)) = (i+1)/n - s/n - y[i]*(1-(i+1)/n);
}
  k1 = k(0,_);
  k2 = k(1,_);
  s=k1[which_max(k1)];
  e=k2[which_max(k2)];
  res= sqrt(n)*std::max(s,e);
return res;
}
