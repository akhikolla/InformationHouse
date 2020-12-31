\name{cont_ks_cdf}
\alias{cont_ks_cdf}
\title{
Computes the cumulative distribution function of the two-sided Kolmogorov-Smirnov statistic when the cdf under the null hypothesis is continuous
}
\description{
Computes the cdf \eqn{P(D_{n} \le q) \equiv P(D_{n} < q)} at a fixed \eqn{q}, \eqn{q\in[0, 1]}, for the one-sample two-sided Kolmogorov-Smirnov statistic, \eqn{D_{n}}, for a given sample size \eqn{n}, when the cdf \eqn{F(x)} under the null hypothesis is continuous.
}
\usage{
cont_ks_cdf(q, n)
}
\arguments{
  \item{q}{
  numeric value between 0 and 1, at which the cdf \eqn{P(D_{n} \le q)} is computed
}
  \item{n}{
  the sample size
  }
}
\details{

Given a random sample \eqn{\{X_{1}, ..., X_{n}\}} of size \code{n} with an empirical cdf \eqn{F_{n}(x)}, the Kolmogorov-Smirnov goodness-of-fit statistic is defined as \eqn{D_{n} = \sup | F_{n}(x) - F(x) | }, where \eqn{F(x)} is the cdf of a prespecified theoretical distribution under the null hypothesis \eqn{H_{0}}, that \eqn{\{X_{1}, ..., X_{n}\}} comes from \eqn{F(x)}.

The function \code{\link{cont_ks_cdf}} implements the FFT-based algorithm proposed by Moscovich and Nadler (2017) to compute the cdf \eqn{P(D_{n} \le q)} at a value \eqn{q}, when \eqn{F(x)} is continuous.
This algorithm ensures a total worst-case run-time of order \eqn{O(n^{2}log(n))} which makes it more efficient and numerically stable than the algorithm proposed by Marsaglia et al. (2003).
The latter is used by many existing packages computing the cdf of \eqn{D_{n}}, e.g., the function \code{\link{ks.test}} in the package \pkg{stats} and the function \code{\link{ks.test}} in the package \pkg{dgof}.
More precisely, in these packages, the exact p-value, \eqn{P(D_{n} \ge q)} is computed only in the case when \eqn{q = d_{n}}, where \eqn{d_{n}} is the value of the KS statistic computed based on a user provided sample \eqn{ \{x_{1}, ..., x_{n} \} }.
Another limitation of the functions \code{\link{ks.test}} is that the sample size should be less than 100, and the computation time is \eqn{O(n^{3})}.
In contrast, the function \code{\link{cont_ks_cdf}} provides results with at least 10 correct digits after the decimal point for sample sizes \eqn{n} up to 100000 and computation time of 16 seconds on a machine with an 2.5GHz Intel Core i5 processor with 4GB RAM, running MacOS X Yosemite.
For \code{n} > 100000, accurate results can still be computed with similar accuracy, but at a higher computation time.
See Dimitrova, Kaishev, Tan (2020), Appendix B for further details and examples.
}
\value{
Numeric value corresponding to \eqn{P(D_{n} \le q)}.
}
\source{
Based on the C++ code available at \url{https://github.com/mosco/crossing-probability} developed by Moscovich and Nadler (2017).
See also Dimitrova, Kaishev, Tan (2020) for more details.
}
\references{
Dimitrina S. Dimitrova, Vladimir K. Kaishev, Senren Tan. (2020) "Computing the Kolmogorov-Smirnov Distribution When the Underlying CDF is Purely Discrete, Mixed or Continuous". Journal of Statistical Software, \bold{95}(10): 1-42. doi:10.18637/jss.v095.i10.

Marsaglia G., Tsang WW., Wang J. (2003). "Evaluating Kolmogorov's Distribution". Journal of Statistical Software, \bold{8}(18), 1-4.

Moscovich A., Nadler B. (2017). "Fast Calculation of Boundary Crossing Probabilities for Poisson Processes". Statistics and Probability Letters, \bold{123}, 177-182.

}


\examples{
## Compute the value for P(D_{100} <= 0.05)

KSgeneral::cont_ks_cdf(0.05, 100)


## Compute P(D_{n} <= q)
## for n = 100, q = 1/500, 2/500, ..., 500/500
## and then plot the corresponding values against q

n<-100
q<-1:500/500
plot(q, sapply(q, function(x) KSgeneral::cont_ks_cdf(x, n)), type='l')

## Compute P(D_{n} <= q) for n = 40, nq^{2} = 0.76 as shown
## in Table 9 of Dimitrova, Kaishev, Tan (2020)

KSgeneral::cont_ks_cdf(sqrt(0.76/40), 40)

}
