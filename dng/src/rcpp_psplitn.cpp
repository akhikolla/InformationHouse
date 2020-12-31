#include <Rcpp.h>
using namespace Rcpp;
//' @describeIn splitn Percentile for the split-normal distribution.
//' @export
// [[Rcpp::export]]
NumericVector psplitn(NumericVector q,NumericVector mu, NumericVector sigma, NumericVector lmd)
{
  int n;
  
  n = q.size();
  mu = rep_len(mu, n);
  sigma = rep_len(sigma, n);
  lmd = rep_len(lmd, n);

  int len;
  len = n;
  NumericVector density(len), out(len);
  NumericVector I0(len),I(len), sign(len);

  for(int a=0;a<len;a++)
  {
    I0[a]=(q[a]<=mu[a]);
    I[a]= 1-I0[a];
    sign[a]=1*I0[a]+lmd[a]*lmd[a]*I[a];
    if(q[a]<=mu[a])
    {
      out[a] =2/(1+lmd[a])*R::pnorm5(q[a],mu[a],sigma[a],1,0);
    }
    else if(q[a]>mu[a])
    {
      out[a] = (1-lmd[a])/(1+lmd[a]) +
        2*lmd[a]/(1+lmd[a])*(R::pnorm5(q[a],mu[a],sigma[a]*lmd[a],1,0)-1/2);
    }
  }

  return out;

}


