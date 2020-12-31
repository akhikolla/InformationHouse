#include <Rcpp.h>
using namespace Rcpp;
//' @describeIn splitn_moments Skewness for the split-normal distribution.
//' @export
// [[Rcpp::export]]
NumericVector splitn_skewness(NumericVector sigma, NumericVector lmd)
{
  IntegerVector a(2);
  int n,j;
  double pi;
  pi = 3.1415926535897932;
  a[0] = sigma.size();
  a[1] = lmd.size();

  if(a[0]==a[1]) {n = a[0];}
  else
  {
    n=a[0];
    if(a[1]>n) n = a[1];
    for(j = a[0];j<n;j++) { sigma[j] = sigma[j-a[0]];}
    for(j = a[1];j<n;j++) { lmd[j] = lmd[j-a[1]];}
  }

  NumericVector skewness(n),nu(n),dm(n);

  for(int i=0;i<n;i++){
    nu[i] = pow(2,0.5)*(lmd[i]-1)*((-4)*pow((lmd[i]-1),2)+pi*(lmd[i]-3)*lmd[i]+pi)*pow(sigma[i],3);
    dm[i] = pow(((-2*pow((lmd[i]-1),2)+pi*lmd[i]*(lmd[i]-1)+pi)*pow(sigma[i],2)),(3/2));

    skewness[i] = nu[i]/dm[i];
  }
  return skewness;
}
