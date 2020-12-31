#include <Rcpp.h>
using namespace Rcpp;
//' @describeIn splitn Randon variables from the split-normal distribution.
//' @export
// [[Rcpp::export]]
NumericVector rsplitn(int n, NumericVector mu, NumericVector sigma, NumericVector lmd)
{
  int i;
  NumericVector u(n),out(n);
  for(i = 0; i<n; i++)
  {
    u[i] = R::runif(0,1);
  }

  mu = rep_len(mu, n);
  sigma = rep_len(sigma, n);
  lmd = rep_len(lmd, n);

  NumericVector mu_long(n),df_long(n),phi_long(n),lmd_long(n);
  NumericVector I0(n),I(n);
  NumericVector p0std(n), y0std(n);

  for(i = 0;i<n;i++)
  {
    I0[i] = (u[i]<=(1/(1+lmd[i])));

    if(I0[i])
    {
      p0std[i] = u[i]*(1+lmd[i])/2;
      y0std[i] = R::qnorm5(p0std[i],mu[i], sigma[i], TRUE, FALSE);
      out[i] = y0std[i]*sigma[i]+mu[i] ;
    }

    else
    {
      p0std[i] = (u[i]-1/(1+lmd[i]))*(1+lmd[i])/(2*lmd[i])+0.5;
      y0std[i] = R::qnorm5(p0std[i], mu[i],sigma[i], TRUE, FALSE);
      out[i] = y0std[i]*(sigma[i]*lmd[i])+mu[i] ;
    }

  }


  return out;
}
