#include <Rcpp.h>
using namespace Rcpp;
//' @describeIn splitn Quantile for the split-normal distribution.
//' @export
// [[Rcpp::export]]
NumericVector qsplitn(NumericVector p,NumericVector mu, NumericVector sigma, NumericVector lmd)
{
  int n;
  n = p.size();
  mu = rep_len(mu, n);
  sigma = rep_len(sigma, n);
  lmd = rep_len(lmd, n);

  NumericVector p0(n),quantile(n);

  for(int i=0;i<n;i++)
  {
    if(p[i]<=(1/(1+lmd[i])))
    {
      p0[i] = (1+lmd[i])*p[i]/2;
      quantile[i] = R::qnorm5(p0[i],mu[i],sigma[i],1,0);

    }

    else
    {
      p0[i] = (p[i]-(1-lmd[i])/(1+lmd[i]))*(1+lmd[i])/(2*lmd[i]);
      quantile[i] = R::qnorm5(p0[i],mu[i],(sigma[i]*lmd[i]),1,0);
    }

  }

  return quantile;

}
