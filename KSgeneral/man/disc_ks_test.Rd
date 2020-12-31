\name{disc_ks_test}
\alias{disc_ks_test}
\title{
Computes the p-value for a one-sample two-sided Kolmogorov-Smirnov test when the cdf under the null hypothesis is purely discrete
}
\description{
Computes the p-value \eqn{P(D_{n} \ge d_{n})}, where \eqn{d_{n}} is the value of the KS test statistic computed based on a data sample \eqn{\{x_{1}, ..., x_{n}\}}, when \eqn{F(x)} is purely discrete, using the Exact-KS-FFT method expressing the p-value as a double-boundary non-crossing probability for a homogeneous Poisson process, which is then efficiently computed using FFT (see Dimitrova, Kaishev, Tan (2020)).
}
\usage{
disc_ks_test(x, y, ..., exact = NULL, tol = 1e-08, sim.size = 1e+06, num.sim = 10)
}
\arguments{
  \item{x}{
  a numeric vector of data sample values \eqn{\{x_{1}, ..., x_{n}\}}.
}
  \item{y}{
   a pre-specified discrete cdf, \eqn{F(x)}, under the null hypothesis.
  Note that \code{y} should be a step function within the class: \code{\link{stepfun}}, of which \code{\link{ecdf}} is a subclass!
}
  \item{\dots}{
    values of the parameters of the cdf, \eqn{F(x)}, specified (as a character string) by \code{y}.
}
  \item{exact}{
  logical variable specifying whether one wants to compute exact p-value \eqn{P(D_{n} \ge d_{n})} using the Exact-KS-FFT method, in which case \code{exact = TRUE} or wants to compute an approximate p-value \eqn{P(D_{n} \ge d_{n})} using the simulation-based algorithm of Wood and Altavela (1978), in which case \code{exact = FALSE}. When \code{exact = NULL} and \code{n <= 100000}, the exact \eqn{P(D_{n} \ge d_{n})} will be computed using the Exact-KS-FFT method. Otherwise, the asymptotic complementary cdf is computed based on Wood and Altavela (1978). By default, \code{exact = NULL}.
}
  \item{tol}{
      the value of \eqn{\epsilon} that is used to compute the values of \eqn{A_{i}} and \eqn{B_{i}}, \eqn{i = 1, ..., n}, as detailed in Step 1 of Section 2.1 in Dimitrova, Kaishev and Tan (2020) (see also (ii) in the Procedure Exact-KS-FFT therein). By default, \code{tol = 1e-08}. Note that a value of \code{NA} or \code{0} will lead to an error!

}
  \item{sim.size}{
   the required number of simulated trajectories in order to produce one Monte Carlo estimate (one MC run) of the asymptotic p-value using the algorithm of Wood and Altavela (1978). By default, \code{sim.size = 1e+06}.

}
  \item{num.sim}{
    the number of MC runs, each producing one estimate (based on \code{sim.size} number of trajectories), which are then averaged in order to produce the final estimate for the asymptotic p-value. This is done in order to reduce the variance of the final estimate. By default, \code{num.sim = 10}.

}
}
\details{

Given a random sample \eqn{\{X_{1}, ..., X_{n}\}} of size \code{n} with an empirical cdf \eqn{F_{n}(x)}, the two-sided Kolmogorov-Smirnov goodness-of-fit statistic is defined as \eqn{D_{n} = \sup | F_{n}(x) - F(x) | }, where \eqn{F(x)} is the cdf of a prespecified theoretical distribution under the null hypothesis \eqn{H_{0}}, that \eqn{\{X_{1}, ..., X_{n}\}} comes from \eqn{F(x)}.

The function \code{\link{disc_ks_test}} implements the Exact-KS-FFT method expressing the p-value as a double-boundary non-crossing probability for a homogeneous Poisson process, which is then efficiently computed using FFT (see Dimitrova, Kaishev, Tan (2020)).
It represents an accurate and fast (run time \eqn{O(n^{2}log(n))}) alternative to the function \code{\link{ks.test}} from the package \pkg{dgof}, which computes a p-value \eqn{P(D_{n} \ge d_{n})}, where \eqn{d_{n}} is the value of the KS test statistic computed based on a user provided data sample \eqn{\{x_{1}, ..., x_{n}\}}, assuming \eqn{F(x)} is purely discrete.

In the function \code{\link{ks.test}}, the p-value for a one-sample two-sided KS test is calculated by combining the approaches of Gleser (1985) and Niederhausen (1981). However, the function \code{\link{ks.test}} due to Arnold and Emerson (2011) only provides exact p-values for \code{n} \eqn{\le} 30, since as noted by the authors, when \code{n} is large, numerical instabilities may occur. In the latter case, \code{\link{ks.test}} uses simulation to approximate p-values, which may be rather slow and inaccurate (see Table 6 of Dimitrova, Kaishev, Tan (2020)).

Thus, making use of the Exact-KS-FFT method, the function \code{\link{disc_ks_test}} provides an exact and highly computationally efficient (alternative) way of computing the p-value \eqn{P(D_{n} \ge d_{n})}, when \eqn{F(x)} is purely discrete.

Lastly, incorporated into the function \code{\link{disc_ks_test}} is the MC simulation-based method of Wood and Altavela (1978) for estimating the asymptotic p-value of \eqn{D_{n}}. The latter method is the default method behind \code{\link{disc_ks_test}} when the sample size \code{n} is \code{n} \eqn{\ge} 100000.


}
\value{
A list with class "htest" containing the following components:
  \item{statistic }{the value of the statistic.}
  \item{p.value }{the p-value of the test.}
  \item{alternative }{"two-sided".}
  \item{data.name }{a character string giving the name of the data.}
}
\references{
Arnold T.A., Emerson J.W. (2011). "Nonparametric Goodness-of-Fit Tests for Discrete Null Distributions". The R Journal, \bold{3}(2), 34-39.

Dimitrina S. Dimitrova, Vladimir K. Kaishev, Senren Tan. (2020) "Computing the Kolmogorov-Smirnov Distribution When the Underlying CDF is Purely Discrete, Mixed or Continuous". Journal of Statistical Software, \bold{95}(10): 1-42. doi:10.18637/jss.v095.i10.

Gleser L.J. (1985). "Exact Power of Goodness-of-Fit Tests of Kolmogorov Type for Discontinuous Distributions". Journal of the American Statistical Association, \bold{80}(392), 954-958.

Niederhausen H. (1981). "Sheffer Polynomials for Computing Exact Kolmogorov-Smirnov and Renyi Type Distributions". The Annals of Statistics, 58-64.

Wood C.L., Altavela M.M. (1978). "Large-Sample Results for Kolmogorov-Smirnov Statistics for Discrete Distributions". Biometrika, \bold{65}(1), 235-239.

}


\seealso{
\code{\link{ks.test}}
}
\examples{
# Comparison of results obtained from dgof::ks.test
# and KSgeneral::disc_ks_test, when F(x) follows the discrete
# Uniform[1, 10] distribution as in Example 3.5 of
# Dimitrova, Kaishev, Tan (2020)

# When the sample size is larger than 100, the
# function dgof::ks.test will be numerically
# unstable

x3 <- sample(1:10, 25, replace = TRUE)
KSgeneral::disc_ks_test(x3, ecdf(1:10), exact = TRUE)
dgof::ks.test(x3, ecdf(1:10), exact = TRUE)
KSgeneral::disc_ks_test(x3, ecdf(1:10), exact = TRUE)$p -
          dgof::ks.test(x3, ecdf(1:10), exact = TRUE)$p

x4 <- sample(1:10, 500, replace = TRUE)
KSgeneral::disc_ks_test(x4, ecdf(1:10), exact = TRUE)
dgof::ks.test(x4, ecdf(1:10), exact = TRUE)
KSgeneral::disc_ks_test(x4, ecdf(1:10), exact = TRUE)$p -
          dgof::ks.test(x4, ecdf(1:10), exact = TRUE)$p

# Using stepfun() to specify the same discrete distribution as defined by ecdf():

steps <- stepfun(1:10, cumsum(c(0, rep(0.1, 10))))
KSgeneral::disc_ks_test(x3, steps, exact = TRUE)

}
