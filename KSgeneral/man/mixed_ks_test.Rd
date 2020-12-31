\name{mixed_ks_test}
\alias{mixed_ks_test}
\title{
Computes the p-value for a one-sample two-sided Kolmogorov-Smirnov test when the cdf under the null hypothesis is mixed
}
\description{
Computes the p-value \eqn{P(D_{n} \ge d_{n})}, where \eqn{d_{n}} is the value of the KS test statistic computed based on a data sample \eqn{\{x_{1}, ..., x_{n}\}}, when \eqn{F(x)} is mixed, using the Exact-KS-FFT method expressing the p-value as a double-boundary non-crossing probability for a homogeneous Poisson process, which is then efficiently computed using FFT (see Dimitrova, Kaishev, Tan (2020)).
}
\usage{
mixed_ks_test(x, jump_points, Mixed_dist, ..., tol = 1e-10)
}

\arguments{
  \item{x}{
    a numeric vector of data sample values \eqn{\{x_{1}, ..., x_{n}\}}.
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

The function \code{\link{mixed_ks_test}} implements the Exact-KS-FFT method expressing the p-value as a double-boundary non-crossing probability for a homogeneous Poisson process, which is then efficiently computed using FFT (see Dimitrova, Kaishev, Tan (2020)).
This algorithm ensures a total worst-case run-time of order \eqn{O(n^{2}log(n))}.

The function \code{\link{mixed_ks_test}} computes the p-value \eqn{P(D_{n} \ge d_{n})}, where \eqn{d_{n}} is the value of the KS test statistic computed based on a user-provided data sample \eqn{\{x_{1}, ..., x_{n}\}}, when \eqn{F(x)} is mixed,

We have not been able to identify alternative, fast and accurate, method (software) that has been developed/implemented when the hypothesized \eqn{F(x)} is mixed.

}
\value{
A list with class "htest" containing the following components:
  \item{statistic }{the value of the statistic.}
  \item{p.value }{the p-value of the test.}
  \item{alternative }{"two-sided".}
  \item{data.name }{a character string giving the name of the data.}
}
\references{
Dimitrina S. Dimitrova, Vladimir K. Kaishev, Senren Tan. (2020) "Computing the Kolmogorov-Smirnov Distribution When the Underlying CDF is Purely Discrete, Mixed or Continuous". Journal of Statistical Software, \bold{95}(10): 1-42. doi:10.18637/jss.v095.i10.

}

\examples{
# Example to compute the p-value of the one-sample two-sided KS test,
# when the underlying distribution is a mixed distribution
# with two jumps at 0 and log(2.5),
# as in Example 3.1 of Dimitrova, Kaishev, Tan (2020)

# Defining the mixed distribution

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
test_data <- c(0,0,0,0,0,0,0.1,0.2,0.3,0.4,
            0.5,0.6,0.7,0.8,log(2.5),log(2.5),
            log(2.5),log(2.5),log(2.5),log(2.5))
KSgeneral::mixed_ks_test(test_data, c(0, log(2.5)),
                         Mixed_cdf_example)


## Compute the p-value of a two-sided K-S test
## when F(x) follows a zero-and-one-inflated
## beta distribution, as in Example 3.3
## of Dimitrova, Kaishev, Tan (2020)

## The data set is the proportion of inhabitants
## living within a 200 kilometer wide costal strip
## in 232 countries in the year 2010

data("Population_Data")
mu <- 0.6189
phi <- 0.6615
a <- mu * phi
b <- (1 - mu) * phi

Mixed_cdf_example <- function(x)
{
     result <- 0
     if (x < 0){
         result <- 0
     }
     else if (x == 0){
         result <- 0.1141
     }
     else if (x < 1){
         result <- 0.1141 + 0.4795 * pbeta(x, a, b)
     }
     else{
         result <- 1
     }

     return (result)
}
KSgeneral::mixed_ks_test(Population_Data, c(0, 1), Mixed_cdf_example)


}
