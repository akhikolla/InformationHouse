#include<Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericVector GoFMI(NumericVector x,double a) { 
  const double pi = 3.14159265358979323846;
  double n = x.size();
  double  e = n/sum(x);    
  NumericVector y = e*x;
  NumericVector res(2);
  NumericVector A1(n*n);
  NumericVector A2(n*n);
  double W1=0;
  double W2=0;
 for (int j=0; j<n;j++){
  for(int k=0;k<n; k++){
    A1[Range(k*n,(k+1)*n-1)]=y-y[k];
    A2[Range(k*n,(k+1)*n-1)]=y+y[k];
  }
 }
 W1=W1+sum(1/(a*a+A1*A1)-1/(a*a+A2*A2)-4*A2/((a*a+A2*A2)*(a*a+A2*A2)));   
 W1=W1+sum((2*a*a-6*A1*A1)/pow((a*a+A1*A1),3)+(2*a*a-6*A2*A2)/pow((a*a+A2*A2),3));
 W2=W2+sum((1+(2*a-A1*A1)/(4*a*a))*exp(-A1*A1/(4*a)));
 W2=W2+sum(((2*a-A2*A2)/(4*a*a)-A2/a-1)*exp(-A2*A2/(4*a)));
 
 W1=a/(2*n)*W1;
 W2=sqrt(pi)/(4*n*sqrt(a))*W2;
 res[0]=W1;
 res[1]=W2;

return res;
}

