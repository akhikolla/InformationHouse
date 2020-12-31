\name{entropy.estimate}
\alias{entropy.estimate}
\title{Vasicek estimate of differential Shannon Entropy}
\description{Computes Vasicek estimate of differential Shannon entropy from a numeric sample.}
\usage{entropy.estimate(x,window)}
\arguments{
  \item{x}{(\code{numeric, vector}) the numeric sample.}
  \item{window}{(\code{numeric, single value}) an integer between 1 and half on the sample size specifying the window size for computing Vasicek estimate. See Details for additional information.}
}
\details{Vasicek estimator of Shannon entropy is defined, for a random sample \eqn{X_1, \dots, X_n}, by
\deqn{\frac{1}{n}\sum_{i=1}^{n} \log (\frac{n}{2m}[X_{(i+m)}-X_{(i-m)}]),}
where \eqn{X_{(i)}} is the order statistic, \eqn{m<(n/2)} is the window size, and \eqn{X_{(i)}=X_{(1)}} for \eqn{i<1} and \eqn{X_{(i)}=X_{(n)}} for \eqn{i>n}. }

\value{ A single numeric value representing the Vasicek estimate of entropy of the sample}
\references{
Vasicek, O., A test for normality based on sample entropy, \emph{Journal of the Royal Statistical Society,} \bold{38(1)}, 54-59 (1976).}
\author{J. Lequesne \email{justine.lequesne@unicaen.fr}}
\seealso{\code{\link{vs.test}} which performs Vasicek-Song goodness-of-fit tests to the specified maximum entropy distribution family.}
\examples{
set.seed(2)
samp <- rnorm(100, mean = 0, s = 1)
entropy.estimate(x = samp, window = 8)
log(2*pi*exp(1))/2 #true value of entropy of normal distribution
}
\keyword{Differential Shannon entropy}
