\name{disc_ks_c_cdf}
\alias{disc_ks_c_cdf}
\title{
Computes the complementary cumulative distribution function of the two-sided Komogorov-Smirnov statistic when the cdf under the null hypothesis is purely discrete
}
\description{

Computes the complementary cdf, \eqn{P(D_{n} \ge q)} at a fixed \eqn{q}, \eqn{q\in[0, 1]}, of the one-sample two-sided Kolmogorov-Smirnov (KS) statistic, when the cdf \eqn{F(x)} under the null hypothesis is purely discrete, using the Exact-KS-FFT method expressing the p-value as a double-boundary non-crossing probability for a homogeneous Poisson process, which is then efficiently computed using FFT (see Dimitrova, Kaishev, Tan (2020)).
Moreover, for comparison purposes, \code{\link{disc_ks_c_cdf}} gives, as an option, the possibility to compute (an approximate value for) the asymptotic \eqn{P(D_{n} \ge q)} using the simulation-based algorithm of Wood and Altavela (1978).
}
\usage{
disc_ks_c_cdf(q, n, y, ..., exact = NULL, tol = 1e-08, sim.size = 1e+06, num.sim = 10)
}
\arguments{
  \item{q}{
     numeric value between 0 and 1, at which the complementary cdf \eqn{P(D_{n}\ge q)} is computed
}
  \item{n}{
    the sample size
}
  \item{y}{
  a pre-specified discrete cdf, \eqn{F(x)} under the null hypothesis.
  Note that \code{y} should be a step function within the class: \code{\link{stepfun}}, of which \code{\link{ecdf}} is a subclass!
}
  \item{\dots}{
    values of the parameters of the cdf, \eqn{F(x)}, specified (as a character string) by \code{y}.
}
  \item{exact}{
  logical variable specifying whether one wants to compute exact p-value \eqn{P(D_{n} \ge q)} using the Exact-KS-FFT method, in which case \code{exact = TRUE} or wants to compute an approximate p-value \eqn{P(D_{n} \ge q)} using the simulation-based algorithm of Wood and Altavela (1978), in which case \code{exact = FALSE}. When \code{exact = NULL} and \code{n <= 100000}, the exact \eqn{P(D_{n} \ge q)} will be computed using the Exact-KS-FFT method. Otherwise, the asymptotic complementary cdf is computed based on Wood and Altavela (1978). By default, \code{exact = NULL}.

}
  \item{tol}{
  the value of \eqn{\epsilon} that is used to compute the values of \eqn{A_{i}} and \eqn{B_{i}}, \eqn{i = 1, ..., n}, as detailed in Step 1 of Section 2.1 in Dimitrova, Kaishev and Tan (2020) (see also (ii) in the Procedure Exact-KS-FFT therein). By default, \code{tol = 1e-08}. Note that a value of \code{NA} or \code{0} will lead to an error!
}
  \item{sim.size}{
  the required number of simulated trajectories in order to produce one Monte Carlo estimate (one MC run) of the asymptotic complementary cdf using the algorithm of Wood and Altavela (1978). By default, \code{sim.size = 1e+06}.

}
  \item{num.sim}{

  the number of MC runs, each producing one estimate (based on \code{sim.size} number of trajectories), which are then averaged in order to produce the final estimate for the asymptotic complementary cdf. This is done in order to reduce the variance of the final estimate. By default, \code{num.sim = 10}.

}
}
\details{

Given a random sample \eqn{\{X_{1}, ..., X_{n}\}} of size \code{n} with an empirical cdf \eqn{F_{n}(x)}, the two-sided Kolmogorov-Smirnov goodness-of-fit statistic is defined as \eqn{D_{n} = \sup | F_{n}(x) - F(x) | }, where \eqn{F(x)} is the cdf of a prespecified theoretical distribution under the null hypothesis \eqn{H_{0}}, that \eqn{\{X_{1}, ..., X_{n}\}} comes from \eqn{F(x)}.

The function \code{\link{disc_ks_c_cdf}} implements the Exact-KS-FFT method, proposed by Dimitrova, Kaishev, Tan (2020) to compute the complementary cdf \eqn{P(D_{n} \ge q)} at a value \eqn{q}, when \eqn{F(x)} is purely discrete.
This algorithm ensures a total worst-case run-time of order \eqn{O(n^{2}log(n))} which makes it more efficient and numerically stable than the only alternative algorithm developed by Arnold and Emerson (2011) and implemented as the function \code{\link{ks.test}} in the package \pkg{dgof}.
The latter only computes a p-value \eqn{P(D_{n} \ge d_{n})}, corresponding to the value of the KS test statistic \eqn{d_{n}} computed based on a user provided sample \eqn{ \{x_{1}, ..., x_{n} \} }.
More precisely, in the package \pkg{dgof} (function \code{\link{ks.test}}), the p-value for a one-sample two-sided KS test is calculated by combining the approaches of Gleser (1985) and Niederhausen (1981). However, the function \code{\link{ks.test}} only provides exact p-values for \code{n} \eqn{\le} 30, since as noted by the authors (see Arnold and Emerson (2011)), when \code{n} is large, numerical instabilities may occur. In the latter case, \code{\link{ks.test}} uses simulation to approximate p-values, which may be rather slow and inaccurate (see Table 6 of Dimitrova, Kaishev, Tan (2020)).

Thus, making use of the Exact-KS-FFT method, the function \code{\link{disc_ks_c_cdf}} provides an exact and highly computationally efficient (alternative) way of computing \eqn{P(D_{n} \ge q)} at a value \eqn{q}, when \eqn{F(x)} is purely discrete.

Lastly, incorporated into the function \code{\link{disc_ks_c_cdf}} is the MC simulation-based method of Wood and Altavela (1978) for estimating the asymptotic complementary cdf of \eqn{D_{n}}. The latter method is the default method behind \code{\link{disc_ks_c_cdf}} when the sample size \code{n} is \code{n} \eqn{\ge} 100000.

}
\value{
Numeric value corresponding to \eqn{P(D_{n} \ge q)}.
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
## Example to compute the exact complementary cdf for D_{n}
## when the underlying cdf F(x) is a binomial(3, 0.5) distribution,
## as shown in Example 3.4 of Dimitrova, Kaishev, Tan (2020)

binom_3 <- stepfun(c(0:3), c(0,pbinom(0:3,3,0.5)))
KSgeneral::disc_ks_c_cdf(0.05, 400, binom_3)

\dontrun{
## Compute P(D_{n} >= q) for n = 100,
## q = 1/5000, 2/5000, ..., 5000/5000, when
## the underlying cdf F(x) is a binomial(3, 0.5) distribution,
## as shown in Example 3.4 of Dimitrova, Kaishev, Tan (2020),
## and then plot the corresponding values against q,
## i.e. plot the resulting complementary cdf of D_{n}

n <- 100
q <- 1:5000/5000
binom_3 <- stepfun(c(0:3), c(0,pbinom(0:3,3,0.5)))
plot(q, sapply(q, function(x) KSgeneral::disc_ks_c_cdf(x, n, binom_3)), type='l')
}

\dontrun{
## Example to compute the asymptotic complementary cdf for D_{n}
## based on Wood and Altavela (1978),
## when the underlying cdf F(x) is a binomial(3, 0.5) distribution,
## as shown in Example 3.4 of Dimitrova, Kaishev, Tan (2020)

binom_3 <- stepfun(c(0: 3), c(0, pbinom(0 : 3, 3, 0.5)))
KSgeneral::disc_ks_c_cdf(0.05, 400, binom_3, exact = FALSE, tol = 1e-08,
sim.size = 1e+06, num.sim = 10)
}

}
