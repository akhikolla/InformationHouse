#include<Rcpp.h>
#include<math.h>
using namespace Rcpp;
// [[Rcpp::export]]
double GoFBH(NumericVector x,double a) { 
  double n = x.size();
  double  e = n/sum(x);    
  NumericVector y=e*x;
  
    double bh=0;
      for (int j=0;j<n; j++){
        
        for(int k=0;k<n; k++){ 
          bh=bh+1/n*((1-y[j])*(1-y[k])/(y[j]+y[k]+a)-(y[j]+y[k])/((y[j]+y[k]+a)*(y[j]+y[k]+a)));
          bh=bh+1/n*((2*y[j]*y[k])/((y[j]+y[k]+a)*(y[j]+y[k]+a)) +(2*y[j]*y[k])/pow((y[j]+y[k]+a),3));      
        } 
      }
    
  return bh;
  }
