#include <Rcpp.h>
using namespace Rcpp;
//' @describeIn splitt Quantile for the split-t distribution.
//' @export
// [[Rcpp::export]]
NumericVector qsplitt(NumericVector p,NumericVector mu, NumericVector df, NumericVector phi, NumericVector lmd)
{
  int n,i;
  n = p.size();
  mu = rep_len(mu, n);
  df = rep_len(df, n);
  phi = rep_len(phi, n);
  lmd = rep_len(lmd, n);

  NumericVector mu_long(n),df_long(n),phi_long(n),lmd_long(n);
  NumericVector I0(n),I(n);
  NumericVector p0std(n), y0std(n);
  NumericVector out(n);


  for(i = 0;i<n;i++)
  {
    I0[i] = (p[i]<=(1/(1+lmd[i])));

    if(I0[i])
    {
      p0std[i] = p[i]*(1+lmd[i])/2;
      y0std[i] = R::qt(p0std[i], df[i], TRUE, FALSE);
      out[i] = y0std[i]*phi[i]+mu[i] ;
    }

    else
    {
      p0std[i] = (p[i]-1/(1+lmd[i]))*(1+lmd[i])/(2*lmd[i])+0.5;
      y0std[i] = R::qt(p0std[i], df[i], TRUE, FALSE);
      out[i] = y0std[i]*(phi[i]*lmd[i])+mu[i] ;
    }

  }
  return out;
}
