#include<Rcpp.h>
#include <math.h>       
using namespace Rcpp;
// [[Rcpp::export]]
NumericVector GoFHM(NumericVector x,double a) { 
  const double pi = 3.14159265358979323846;

  double n = x.size();
  double  e = n/sum(x);    
  NumericVector y=e*x;
  NumericVector res(2);
  NumericMatrix Ajkm(n,n);
  NumericMatrix Ajkp(n,n);
  double a2=a*a;
  double T1=0;
  double T2=0;
  for (int j = 0; j<n; j++){
   for(int k=0; k<n; k++){
     Ajkm(_,k)=y-y[k];
     Ajkp(_,k)=y+y[k];
   }
 }
 T1=(a/n)*(sum(1/(a2+pow(Ajkm,2)))+sum(1/(a2+pow(Ajkp,2))));
 T2=1/(2*n)*sqrt(pi/a)*(sum(exp(-pow(Ajkm,2)/(4*a))+exp(-pow(Ajkp,2)/(4*a))));
//Rcout<<T1<<" "<<T2<<std::endl;
 for(int l=0;l<n; l++){

   T1=T1-2*a/pow(n,2)*(sum(1/(a2+pow(Ajkm-y[l],2))+1/(a2+pow(Ajkm+y[l],2))));
   T2=T2-1/pow(n,2)*sqrt(pi/a)*sum(exp(-pow(Ajkm-y[l],2)/(4*a))+exp(-pow(Ajkm+y[l],2)/(4*a)));
   for (int m=0; m<n;m++){
     T1=T1+a/pow(n,3)*sum(1/(pow(a,2)+(Ajkm-Ajkm(l,m))*(Ajkm-Ajkm(l,m))));
     T1=T1+a/pow(n,3)*sum(1/(pow(a,2)+(Ajkm+Ajkm(l,m))*(Ajkm+Ajkm(l,m))));
     T2=T2+1/(2*pow(n,3))*sqrt(pi/a)*sum(exp(-(Ajkm-Ajkm(l,m))*(Ajkm-Ajkm(l,m))/(4*a))+exp(-(Ajkm+Ajkm(l,m))*(Ajkm+Ajkm(l,m))/(4*a)));

  }
}
res[0]=T1;
res[1]=T2;
return res;
}
