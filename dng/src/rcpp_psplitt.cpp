#include <Rcpp.h>
using namespace Rcpp;
//' @describeIn splitt Percentile for the split-t distribution.
//' @export
// [[Rcpp::export]]
NumericVector psplitt(NumericVector q, NumericVector mu, NumericVector df, NumericVector phi, NumericVector lmd)
{
  double ibeta0;
  double pbeta0;
  int n,i;

  n = q.size();
  mu = rep_len(mu, n);
  df = rep_len(df, n);
  phi = rep_len(phi, n);
  lmd = rep_len(lmd, n);
  
  NumericVector I0(n),I(n), sign(n), sign2(n);
  NumericVector A(n),BetaRegUpper(n);
  NumericVector out(n);

  for(i = 0;i<n;i++){
    I0[i] = (q[i]<=mu[i]);
    I[i]  = 1-I0[i];
    sign[i]  = 1*I0[i]+lmd[i]*I[i];
    sign2[i] = -1*I0[i]+1*I[i];

    A[i] = df[i]*pow(sign[i],2)*pow(phi[i],2)/(df[i]*pow(sign[i],2)*pow(phi[i],2)+pow((q[i]-mu[i]),2));

    pbeta0 = ::Rf_pbeta(A[i],df[i]*0.5, 0.5, TRUE, TRUE);
    ibeta0 = exp(pbeta0);

    BetaRegUpper[i] = 1-ibeta0;

    out[i] = (1/(1+lmd[i]) + sign[i]*sign2[i]/(1+lmd[i])*BetaRegUpper[i]);

  }

  return out;
}
