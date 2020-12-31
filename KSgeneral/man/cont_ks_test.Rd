\name{cont_ks_test}
\alias{cont_ks_test}
\title{
Computes the p-value for a one-sample two-sided Kolmogorov-Smirnov test when the cdf under the null hypothesis is continuous
}
\description{
Computes the p-value \eqn{P(D_{n} \ge d_{n}) \equiv P(D_{n} > d_{n})}, where \eqn{d_{n}} is the value of the KS test statistic computed based on a data sample \eqn{\{x_{1}, ..., x_{n}\}}, when \eqn{F(x)} is continuous.

}
\usage{
cont_ks_test(x, y, ...)
}
\arguments{
  \item{x}{
  a numeric vector of data sample values \eqn{\{x_{1}, ..., x_{n}\}}.
}
  \item{y}{
  a pre-specified continuous cdf, \eqn{F(x)} under the null hypothesis. Note that \code{y} should be a character string naming a continuous cumulative distribution function such as \code{\link{pexp}}, \code{\link{pnorm}}, etc.
  Only continuous cdfs are valid!
}
  \item{\dots}{
  values of the parameters of the cdf, \eqn{F(x)} specified (as a character string) by \code{y}.
}
}
\details{
Given a random sample \eqn{\{X_{1}, ..., X_{n}\}} of size \code{n} with an empirical cdf \eqn{F_{n}(x)}, the two-sided Kolmogorov-Smirnov goodness-of-fit statistic is defined as \eqn{D_{n} = \sup | F_{n}(x) - F(x) | }, where \eqn{F(x)} is the cdf of a prespecified theoretical distribution under the null hypothesis \eqn{H_{0}}, that \eqn{\{X_{1}, ..., X_{n}\}} comes from \eqn{F(x)}.

The function \code{\link{cont_ks_test}} implements the FFT-based algorithm proposed by Moscovich and Nadler (2017) to compute the p-value \eqn{P(D_{n} \ge d_{n})}, where \eqn{d_{n}} is the value of the KS test statistic computed based on a user provided data sample \eqn{\{x_{1}, ..., x_{n}\}}, assuming \eqn{F(x)} is continuous.
This algorithm ensures a total worst-case run-time of order \eqn{O(n^{2}log(n))} which makes it more efficient and numerically stable than the algorithm proposed by Marsaglia et al. (2003).
The latter is used by many existing packages computing the cdf of \eqn{D_{n}}, e.g., the function \code{\link{ks.test}} in the package \pkg{stats} and the function \code{\link{ks.test}} in the package \pkg{dgof}.
A limitation of the functions \code{\link{ks.test}} is that the sample size should be less than 100, and the computation time is \eqn{O(n^{3})}.
In contrast, the function \code{\link{cont_ks_test}} provides results with at least 10 correct digits after the decimal point for sample sizes \eqn{n} up to 100000 and computation time of 16 seconds on a machine with an 2.5GHz Intel Core i5 processor with 4GB RAM, running MacOS X Yosemite.
For \code{n} > 100000, accurate results can still be computed with similar accuracy, but at a higher computation time.
See Dimitrova, Kaishev, Tan (2020), Appendix C for further details and examples.

}
\value{
A list with class "htest" containing the following components:

  \item{statistic }{the value of the statistic.}
  \item{p.value }{the p-value of the test.}
  \item{alternative }{"two-sided".}
  \item{data.name }{a character string giving the name of the data.}

}
\source{
Based on the C++ code available at \url{https://github.com/mosco/crossing-probability} developed by Moscovich and Nadler (2017).
See also Dimitrova, Kaishev, Tan (2020) for more details.
}

\references{
Dimitrina S. Dimitrova, Vladimir K. Kaishev, Senren Tan. (2020) "Computing the Kolmogorov-Smirnov Distribution When the Underlying CDF is Purely Discrete, Mixed or Continuous". Journal of Statistical Software, \bold{95}(10): 1-42. doi:10.18637/jss.v095.i10.

Moscovich A., Nadler B. (2017). "Fast Calculation of Boundary Crossing Probabilities for Poisson Processes". Statistics and Probability Letters, \bold{123}, 177-182.
}

\examples{

## Comparing the p-values obtained by stat::ks.test
## and KSgeneral::cont_ks_test

x<-abs(rnorm(100))
p.kt <- ks.test(x, "pexp", exact = TRUE)$p
p.kt_fft <- KSgeneral::cont_ks_test(x, "pexp")$p
abs(p.kt-p.kt_fft)


}
