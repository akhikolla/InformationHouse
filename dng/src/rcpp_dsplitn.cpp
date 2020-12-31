#include <Rcpp.h>
using namespace Rcpp;
//' Split-normal distribution
//'
//' Density distribution function, quantile function and random generation function for
//' the split normal distribution.
//'
//' The random ' variable y follows a split-normal distribution, y~N(\eqn{\mu}, '
//' \eqn{\sigma}, \eqn{\lambda}), which has density: \deqn{1/(1+\lambda)\sigma ' \sqrt(2/\pi)
//' exp{-(y-\mu)*2/2\sigma^2}, if y<=\mu} ' \deqn{1/(1+\lambda)\sigma \sqrt(2/\pi)
//' exp{-(y-\mu)*2/2\sigma^2 \lambda^2}, ' if y>\mu} where \eqn{\sigma>0} and
//' \eqn{\lambda>0}. The Split-normal ' distribution reduce to normal distribution when
//' \eqn{\lambda=1}.
//'
//' @name splitn
//' @param x vector of quantiles.
//' @param mu vector of location parameter. (The mode of the density)
//' @param sigma vector of standard deviations.
//' @param lmd vector of skewness parameters (>0). If is 1, reduced to
//' symmetric normal distribution.
//' @param p vector of probability.
//' @param q vector of quantiles.
//' @param n number of observations. If length(n) > 1, the length is taken to be the number required.
//' @param logarithm logical; if TRUE, probabilities p are given as log(p).
//' @return \code{dsplitn} gives the density; \code{psplitn} gives the percentile;
//' \code{qsplitn} gives the quantile; and \code{rsplitn} gives the random
//' variables. Invalid arguments will result in return value NaN, with a warning.
//'
//' The numerical arguments other than n are recycled to the length of the
//' result. Only the first elements of the logical arguments are used.
//'
//' @author Feng Li, Jiayue Zeng
//' @seealso \code{\link{splitn_mean}()},
//' \code{\link{splitn_var}()},\code{\link{splitn_skewness}()} and
//' \code{\link{splitn_kurtosis}()} for numerical characteristics of the
//' split-normal distribution.
//'
//' @references
//' Villani, M., & Larsson, R. (2006) The Multivariate Split Normal
//' Distribution and Asymmetric Principal Components Analysis. Sveriges
//' Riksbank Working Paper Series, No. 175.
//'
//' @keywords distribution asymmetric normal
//' @examples
//'
//' n <- 3
//' mu <- c(0,1,2)
//' sigma <- c(1,2,3)
//' lmd <- c(1,2,3)
//'
//' q0 <- rsplitn(n, mu, sigma, lmd)
//' d0 <- dsplitn(q0, mu, sigma, lmd, logarithm = FALSE)
//' p0 <- psplitn(q0, mu, sigma, lmd)
//' q1 <- qsplitn(p0,mu, sigma, lmd)
//' all.equal(q0, q1)
//' @export
// [[Rcpp::export]]
NumericVector dsplitn(NumericVector x, NumericVector mu, NumericVector sigma, NumericVector lmd, bool logarithm)
{
  int n;

  n = x.size();
  mu = rep_len(mu, n);
  sigma = rep_len(sigma, n);
  lmd = rep_len(lmd, n);

  int len;
  double pi;
  pi = 3.1415926535897932;
  len = n;
  NumericVector densitq(len),out(len);
  NumericVector I0(len),I(len), sign(len);

  for(int a=0;a<len;a++)
  {
    I0[a]=(x[a]<=mu[a]);
    I[a]= 1-I0[a];
    sign[a]=1*I0[a]+lmd[a]*lmd[a]*I[a];
    densitq[a] = sqrt(2/pi)*
      exp(-pow((x[a]-mu[a]),2)/(2*sigma[a]*sigma[a]*sign[a]))/
        ((1+lmd[a])*sigma[a]);
  }

  if(logarithm)
  {
    for(int i = 0;i<len;i++)
    { out[i] = exp(densitq[i]);   }
  }
  else {out = densitq;}
  return out;
}
