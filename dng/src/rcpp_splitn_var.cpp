#include <Rcpp.h>
using namespace Rcpp;
//' @describeIn splitn_moments Variance for the split-normal distribution.
//' @export
// [[Rcpp::export]]
NumericVector splitn_var(NumericVector sigma, NumericVector lmd)
{
  IntegerVector a(2);
  int n,i,j;
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

  NumericVector var(n),b(n);

  for(i=0;i<n;i++){
    b[i] = -2/pi*pow((lmd[i]-1),2)+lmd[i]*(lmd[i]-1)+1;
    var[i] = b[i]*pow(sigma[i],2);
  }
  return var;
}
