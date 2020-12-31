#include <Rcpp.h>
using namespace Rcpp;
//' @describeIn splitn_moments Kurtosis for the split-normal distribution.
//' @export
// [[Rcpp::export]]
NumericVector splitn_kurtosis(NumericVector lmd)
{
  int n;
  double pi;
  n =lmd.size();
  pi = 3.1415926535897932;

  NumericVector kurtosis(n);
  NumericVector k1(n),k2(n),k3(n);

  for(int i=0;i<n;i++){
    k1[i] = pow((lmd[i]-1),2);
    k2[i] = 8*(pi-3)*lmd[i]*lmd[i]+3*pow((pi-4),2)+8*(pi-3);
    k3[i] = pow((-2*pow((lmd[i]-1),2)+pi*lmd[i]*(lmd[i]-1)+pi),2);
    kurtosis[i] = k1[i]*k2[i]/k3[i];
  }
  return kurtosis;
}
