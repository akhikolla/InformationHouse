#include <Rcpp.h>
using namespace Rcpp;
//' Moments of the split normal distribution
//'
//' Computing the mean, variance, skewness and kurtosis for the split-normal
//' distribution.
//'
//'
//' @name splitn_moments
//' @param mu vector of location parameter. (The mode of the density)
//' @param sigma vector of standard deviations.
//' @param lmd vector of skewness parameters (>0). If is 1, reduce to normal
//' distribution.
//' @return \code{splitn_mean} gives the mean.  \code{splitn_var} gives the
//' variance.  \code{splitn_skewness} gives the skewness.
//' \code{splitn_kurtosis} gives the kurtosis.  (\code{splitn_mean},
//' \code{splitn_var},\code{splitn_skeness} and \code{splitn_kurtosis} are all
//' vectors.
//'
//' @author Feng Li, Jiayue Zeng
//'
//' @seealso \code{\link{psplitn}()} \code{\link{dsplitn}()} \code{\link{qsplitn}()} and
//' \code{\link{rsplitn}()} for the split-normal distribution.
//'
//' @references
//' Villani, M., & Larsson, R. (2006) The Multivariate Split Normal
//' Distribution and Asymmetric Principal Components Analysis. Sveriges
//' Riksbank Working Paper Series, No. 175.
//'
//' @keywords distribution asymmetric normal
//' @examples
//'
//' mu <- c(0,1,2)
//' sigma <- c(0.5,1,2)
//' lmd <- c(1,2,3)
//'
//' mean0 <- splitn_mean(mu, sigma, lmd)
//' var0 <- splitn_var(sigma, lmd)
//' skewness0 <- splitn_skewness(sigma, lmd)
//' kurtosis0 <- splitn_kurtosis(lmd)
//' @export
// [[Rcpp::export]]
NumericVector splitn_mean(NumericVector mu, NumericVector sigma, NumericVector lmd)
{
  IntegerVector a(3);
  int n,i,j;
  double pi;
  pi = 3.1415926535897932;
  a[0] = mu.size();
  a[1] = sigma.size();
  a[2] = lmd.size();

  if(a[0]==a[1] && a[0]==a[2]) {n = a[0];}
  else
  {
    n=a[0];
    for(i = 1;i<=2;i++)   { if(a[i]>n) n = a[i];}
    for(j = a[0];j<n;j++) { mu[j] = mu[j-a[0]];}
    for(j = a[1];j<n;j++) { sigma[j] = sigma[j-a[1]];}
    for(j = a[2];j<n;j++) { lmd[j] = lmd[j-a[2]];}
  }

  NumericVector mean(n);

  for(int i=0;i<n;i++){
    mean[i] = sqrt(2/pi)*(lmd[i]-1)*sigma[i]+mu[i];
  }
  return mean;
}
