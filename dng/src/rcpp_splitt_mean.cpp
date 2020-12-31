#include <Rcpp.h>
using namespace Rcpp;
//' Moments of the split-t distribution
//'
//' Computing the mean, variance, skewness and kurtosis for the split student-t
//' distribution.
//'
//'
//' @name splitt_moments
//' @param mu vector of location parameter. (The mode of the density)
//' @param df degrees of freedom (> 0, can be non-integer). df = Inf is allowed.
//' @param phi vector of scale parameters (> 0).
//' @param lmd vector of skewness parameters (> 0). If is 1, reduced to
//' symmetric student t distribution.
//' @return \code{splitt_mean} gives the mean. \code{splitt_var} gives the
//' variance. \code{splitt_skewness} gives the skewness. \code{splitt_kurtosis}
//' gives the kurtosis. (\code{splitt_mean},
//' \code{splitt_var},\code{splitt_skeness} and \code{splitt_kurtosis} are all
//' vectors.)
//'
//' Invalid arguments will result in return value NaN, with a warning.
//' @author Feng Li, Jiayue Zeng
//' @seealso \code{\link{dsplitt}()}, \code{\link{psplitt}()},
//' \code{\link{qsplitt}()} and \code{\link{rsplitt}()} for the split-t
//' distribution.
//'
//' @references
//' Li, F., Villani, M., & Kohn, R. (2010). Flexible modeling of
//' conditional distributions using smooth mixtures of asymmetric student t
//' densities. Journal of Statistical Planning & Inference, 140(12), 3638-3654.
//' @keywords distribution asymmetric student-t
//' @examples
//'
//' mu <- c(0,1,2)
//' df <- rep(10,3)
//' phi <- c(0.5,1,2)
//' lmd <- c(1,2,3)
//'
//' mean0 <- splitt_mean(mu, df, phi, lmd)
//' var0 <- splitt_var(df, phi, lmd)
//' skewness0 <- splitt_skewness(df, phi, lmd)
//' kurtosis0 <- splitt_kurtosis(df, phi, lmd)
//' @export
// [[Rcpp::export]]
NumericVector splitt_mean(NumericVector mu, NumericVector df, NumericVector phi, NumericVector lmd)
{
  IntegerVector a(4);
  int n,i,j;
  a[0] = mu.size();
  a[1] = df.size();
  a[2] = phi.size();
  a[3] = lmd.size();

  if(a[0]==a[1] && a[0]==a[2] && a[0]==a[3]) {n = a[0];}
  else
  {
    n=a[0];
    for(i = 1;i<=3;i++)   { if(a[i]>n) n = a[i];}
    for(j = a[0];j<n;j++) { mu[j] = mu[j-a[0]];}
    for(j = a[1];j<n;j++) { df[j] = df[j-a[1]];}
    for(j = a[2];j<n;j++) { phi[j] = phi[j-a[2]];}
    for(j = a[3];j<n;j++) { lmd[j] = lmd[j-a[3]];}
  }

  NumericVector h(n),mean(n);
  NumericVector beta0(n);

  for(int i=0;i<n;i++){
    beta0[i]=R::beta(df[i]*0.5,0.5);
    h[i] = 2*pow(df[i],0.5)*phi[i]*(lmd[i]-1)/((df[i]-1)*beta0[i]);
    mean[i] = mu[i]+h[i];
  }
  return mean;
}
