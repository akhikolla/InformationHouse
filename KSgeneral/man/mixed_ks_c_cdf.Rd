\name{mixed_ks_c_cdf}
\alias{mixed_ks_c_cdf}
\title{
Computes the complementary cumulative distribution function of the two-sided Kolmogorov-Smirnov statistic when the cdf under the null hypothesis is mixed
}
\description{
Computes the complementary cdf, \eqn{P(D_{n} \ge q)} at a fixed \eqn{q}, \eqn{q\in[0, 1]}, of the one-sample two-sided Kolmogorov-Smirnov statistic,  when the cdf \eqn{F(x)} under the null hypothesis is mixed, using the Exact-KS-FFT method expressing the p-value as a double-boundary non-crossing probability for a homogeneous Poisson process, which is then efficiently computed using FFT (see Dimitrova, Kaishev, Tan (2020)).
}
\usage{
mixed_ks_c_cdf(q, n, jump_points, Mixed_dist, ..., tol = 1e-10)
}
\arguments{
  \item{q}{
  numeric value between 0 and 1, at which the complementary cdf \eqn{P(D_{n} \ge q)} is computed
}
  \item{n}{
      the sample size
}
  \item{jump_points}{
  a numeric vector containing the points of (jump) discontinuity, i.e. where the underlying cdf \eqn{F(x)} has jump(s)
}
  \item{Mixed_dist}{
    a pre-specified (user-defined) mixed cdf, \eqn{F(x)}, under the null hypothesis.
}
  \item{\dots}{
      values of the parameters of the cdf, \eqn{F(x)} specified (as a character string) by \code{Mixed_dist}.
}
  \item{tol}{

   the value of \eqn{\epsilon} that is used to compute the values of \eqn{A_{i}} and \eqn{B_{i}}, \eqn{i = 1, ..., n}, as detailed in Step 1 of Section 2.1 in Dimitrova, Kaishev and Tan (2020) (see also (ii) in the Procedure Exact-KS-FFT therein). By default, \code{tol = 1e-10}. Note that a value of \code{NA} or \code{0} will lead to an error!

}
}
\details{
Given a random sample \eqn{\{X_{1}, ..., X_{n}\}} of size \code{n} with an empirical cdf \eqn{F_{n}(x)}, the Kolmogorov-Smirnov goodness-of-fit statistic is defined as \eqn{D_{n} = \sup | F_{n}(x) - F(x) | }, where \eqn{F(x)} is the cdf of a prespecified theoretical distribution under the null hypothesis \eqn{H_{0}}, that \eqn{\{X_{1}, ..., X_{n}\}} comes from \eqn{F(x)}.

The function \code{\link{mixed_ks_c_cdf}} implements the Exact-KS-FFT method, proposed by Dimitrova, Kaishev, Tan (2020) to compute the complementary cdf \eqn{P(D_{n} \ge q)} at a value \eqn{q}, when \eqn{F(x)} is mixed.
This algorithm ensures a total worst-case run-time of order \eqn{O(n^{2}log(n))}.

We have not been able to identify alternative, fast and accurate, method (software) that has been developed/implemented when the hypothesized \eqn{F(x)} is mixed.

}
\value{
Numeric value corresponding to \eqn{P(D_{n} \ge q)}.
}
\references{
Dimitrina S. Dimitrova, Vladimir K. Kaishev, Senren Tan. (2020) "Computing the Kolmogorov-Smirnov Distribution When the Underlying CDF is Purely Discrete, Mixed or Continuous". Journal of Statistical Software, \bold{95}(10): 1-42. doi:10.18637/jss.v095.i10.

}


\examples{
# Compute the complementary cdf of D_{n}
# when the underlying distribution is a mixed distribution
# with two jumps at 0 and log(2.5),
# as in Example 3.1 of Dimitrova, Kaishev, Tan (2020)

## Defining the mixed distribution

Mixed_cdf_example <- function(x)
{
     result <- 0
     if (x < 0){
         result <- 0
     }
     else if (x == 0){
         result <- 0.5
     }
     else if (x < log(2.5)){
         result <- 1 - 0.5 * exp(-x)
     }
     else{
         result <- 1
     }

     return (result)
 }

KSgeneral::mixed_ks_c_cdf(0.1, 25, c(0, log(2.5)), Mixed_cdf_example)


\dontrun{
## Compute P(D_{n} >= q) for n = 5,
## q = 1/5000, 2/5000, ..., 5000/5000
## when the underlying distribution is a mixed distribution
## with four jumps at 0, 0.2, 0.8, 1.0,
## as in Example 2.8 of Dimitrova, Kaishev, Tan (2020)

n <- 5
q <- 1:5000/5000

Mixed_cdf_example <- function(x)
{
  result <- 0
  if (x < 0){
    result <- 0
  }
  else if (x == 0){
    result <- 0.2
  }
  else if (x < 0.2){
    result <- 0.2 + x
  }
  else if (x < 0.8){
    result <- 0.5
  }
  else if (x < 1){
    result <- x - 0.1
  }
  else{
    result <- 1
  }

  return (result)
}

plot(q, sapply(q, function(x) KSgeneral::mixed_ks_c_cdf(x, n,
     c(0, 0.2, 0.8, 1.0), Mixed_cdf_example)), type='l')

}

}
