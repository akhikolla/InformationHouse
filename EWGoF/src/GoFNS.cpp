#include <Rcpp.h>
#include<math.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericVector GoFNS(double t, double n,double m) {
  NumericVector res(m);
  double p_r,q_r,Q_r,d_Q_r = 0;
  double d2_Q_r, d3_Q_r,d4_Q_r=0;
  for(int i=0 ; i<m; ++i){
   p_r = (t+i)/(n+1);
   q_r = 1-p_r;
   Q_r = log(-log(1-p_r));
   d_Q_r = -1/((1-p_r)*log(1-p_r));
   d2_Q_r =  d_Q_r/q_r - d_Q_r*d_Q_r;
   d3_Q_r = d2_Q_r/q_r + d_Q_r/(q_r*q_r) - 2*(d_Q_r*d2_Q_r);
   d4_Q_r = d3_Q_r/q_r  + 2*d2_Q_r/(q_r*q_r);
   d4_Q_r += 2*d_Q_r/(q_r*q_r*q_r)- 2*(d2_Q_r*d2_Q_r+d_Q_r*d3_Q_r);
   res[i] = Q_r + (p_r*q_r/(2*(n+2)))*(d2_Q_r);
   res[i] += q_r*p_r/((n+2)*(n+2))*(1/3*(q_r - p_r)*d3_Q_r + 1/8*p_r*q_r*d4_Q_r);
}
NumericVector temp=Rcpp::wrap(res);
return res;}
