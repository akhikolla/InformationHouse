#include <Rcpp.h>
using namespace Rcpp;
//' Split-t distribution
//'
//' Density, distribution function, quantile function and random generation for the normal
//' distribution for the split student-t distribution.
//'
//' The random variable y follows a split-t distribution with \eqn{\nu}>0
//' degrees of freedom, y~t(\eqn{\mu}, \eqn{\phi}, \eqn{\lambda}, \eqn{\nu}),
//' if its density function is of the form
//'
//' \deqn{C K(\mu, \phi, \nu,)I(y\leq\mu) + C K(\mu, \lambda \phi,
//' \nu)I(y>\mu), } where, \deqn{K(\mu, \phi, \nu,) =[\nu/(\nu+(y-\mu)^2 /\phi
//' ^2)]^{(\nu+1)/2} } is the kernel of a student \eqn{t} density with variance
//' \eqn{\phi ^2\nu/(\nu-2)} and \deqn{c = 2[(1+\lambda)\phi (\sqrt \nu)
//' Beta(\nu/2,1/2)]^{-1} }is the normalization constant.
//'
//' @name splitt
//'
//' @param x vector of quantiles.
//' @param mu vector of location parameter. (The mode of the density)
//' @param df degrees of freedom (> 0, can be non-integer). df = Inf is also allowed.
//' @param phi vector of scale parameters (>0).
//' @param lmd vector of skewness parameters (>0). If is 1, reduced to the
//' symmetric student t distribution.
//' @param p vector of probability.
//' @param q vector of quantiles.
//' @param n number of observations. If length(n) > 1, the length is taken to be the number required.
//' @param logarithm logical; if TRUE, probabilities p are given as log(p).
//' @return \code{dsplitt} gives the density; \code{psplitt} gives the percentile;
//' \code{qsplitt} gives the quantile; and \code{rsplitt} gives the random
//' variables. Invalid arguments will result in return value NaN, with a warning.
//'
//' The numerical arguments other than n are recycled to the length of the
//' result. Only the first elements of the logical arguments are used.
//'
//' @author Feng Li, Jiayue Zeng
//' @seealso \code{\link{splitt_mean}()},
//' \code{\link{splitt_var}()},\code{\link{splitt_skewness}()} and
//' \code{\link{splitt_kurtosis}()} for numerical characteristics of the
//' Split-t distribution.
//'
//' @references
//' Li, F., Villani, M., & Kohn, R. (2010). Flexible modeling of
//' conditional distributions using smooth mixtures of asymmetric student t
//' densities. Journal of Statistical Planning & Inference, 140(12), 3638-3654.
//' @keywords distribution asymmetric student-t
//'
//' @examples
//'
//' n <- 3
//' mu <- c(0,1,2)
//' df <- rep(10,3)
//' phi <- c(0.5,1,2)
//' lmd <- c(1,2,3)
//'
//' q0 <- rsplitt(n, mu, df, phi, lmd)
//' d0 <- dsplitt(q0, mu, df, phi, lmd, logarithm = FALSE)
//' p0 <- psplitt(q0, mu, df, phi, lmd)
//' q1 <- qsplitt(p0,mu, df, phi, lmd)
//' all.equal(q0, q1)
//'
//' @export
// [[Rcpp::export]]
NumericVector dsplitt(NumericVector x,NumericVector mu, NumericVector df, NumericVector phi, NumericVector lmd, bool logarithm)
{
  int n,i;

  n = x.size();
  mu = rep_len(mu, n);
  df = rep_len(df, n);
  phi = rep_len(phi, n);
  lmd = rep_len(lmd, n);

  NumericVector I0(n),I(n);
  NumericVector sign(n),densitylog(n),out(n);
  NumericVector lbeta0(n);

    for(i = 0;i<n;i++)
    {
      lbeta0[i] = ::Rf_lbeta(0.5*df[i],0.5);
      I0[i] = (x[i]<=mu[i]); // Logical values. 1, if y <= mu; 0, if y >mu.
      I[i] = (x[i]>mu[i]); //Logical values. 1, if y > mu; 0, if y <= mu.
      sign[i] = 1*I0[i]+lmd[i]*I[i]; // sign = 1 if y<=mu; sign = lmd.^2 if y>2
      densitylog[i] = (std::log(2)+(1+df[i])/2*(std::log(df[i])-std::log(df[i]+pow((-mu[i]+x[i]),2)/(pow(phi[i],2)*pow(sign[i],2))))-std::log(phi[i])-std::log(df[i])/2-lbeta0[i]-std::log(1+lmd[i]));
      out[i] = densitylog[i];
    }


  if(!logarithm)
  {
    for(i = 0;i<n;i++){ out[i] = exp(densitylog[i]);   }
  }

  return out;
}
