#include <Rcpp.h>
using namespace Rcpp;
//' @describeIn splitt_moments Skewness for the split-t distribution.
//' @export
// [[Rcpp::export]]
NumericVector splitt_skewness(NumericVector df, NumericVector phi, NumericVector lmd)
{
  IntegerVector a(3);
  int n,i,j;
  a[0] = df.size();
  a[1] = phi.size();
  a[2] = lmd.size();

  if(a[0]==a[1] && a[0]==a[2]) {n=a[0];}
  else
  {
    n=a[0];
    for(i=1;i<=2;i++)   { if(a[i]>n) n=a[i];}
    for(j=a[0];j<n;j++) { df[j]=df[j-a[0]];}
    for(j=a[1];j<n;j++) { phi[j]=phi[j-a[1]];}
    for(j=a[2];j<n;j++) { lmd[j]=lmd[j-a[2]];}
  }

  NumericVector h(n),var(n),m3(n),skewness(n);
  NumericVector beta0(n);

  for(int i=0;i<n;i++){
    beta0[i] = R::beta(df[i]*0.5,0.5);
    h[i]     = 2*pow(df[i],0.5)*phi[i]*(lmd[i]-1)/((df[i]-1)*beta0[i]);
    var[i]   = (1+pow(lmd[i],3))/(1+lmd[i])*df[i]/(df[i]-2)*phi[i]*phi[i]-h[i]*h[i];
    m3[i]    = (2*pow(h[i],3)+2*h[i]*pow(phi[i],2)*(pow(lmd[i],2)+1)*df[i]/(df[i]-3)
             -3*h[i]*pow(phi[i],2)*(pow(lmd[i],3)+1)/(lmd[i]+1)*df[i]/(df[i]-2));

    skewness[i] = m3[i]/pow(var[i],1.5);
  }
  return skewness;
}
