\name{KSgeneral-package}
\alias{KSgeneral-package}
\docType{package}
\title{
\packageTitle{KSgeneral}
}
\description{
The one-sample two-sided Kolmogorov-Smirnov (KS) statistic is one of the most popular goodness-of-fit test statistics that is used to measure how well the distribution of a random sample agrees with a prespecified theoretical distribution.
Given a random sample \eqn{\{X_{1}, ..., X_{n}\}} of size \eqn{n} with an empirical cdf \eqn{F_{n}(x)}, the two-sided KS statistic is defined as
\eqn{D_{n} = \sup | F_{n}(x) - F(x) | }, where \eqn{F(x)} is the cdf of the prespecified theoretical distribution under the null hypothesis \eqn{H_{0}}, that \eqn{ \{ X_{1}, ..., X_{n} \} } comes from \eqn{F(x)}.



The package \pkg{KSgeneral} implements a novel, accurate and efficient Fast Fourier Transform (FFT)-based method, referred as Exact-KS-FFT method to compute the complementary cdf,
\eqn{P(D_{n} \ge q)}, at a fixed \eqn{q\in[0, 1]} for a given (hypothezied) purely discrete, mixed or continuous underlying cdf \eqn{F(x)}, and arbitrary, possibly large sample size \eqn{n}.
A plot of the complementary cdf \eqn{P(D_{n} \ge q)}, \eqn{0 \le q \le 1}, can also be produced.

In other words, the package computes the p-value, \eqn{P(D_{n} \ge q)} for any fixed critical level \eqn{q\in[0, 1]}.
If a data sample, \eqn{\{x_{1}, ..., x_{n}\}} is supplied, \pkg{KSgeneral} computes the p-value \eqn{P(D_{n} \ge d_{n})}, where \eqn{d_{n}} is the value of the KS test statistic computed based on \eqn{\{x_{1}, ..., x_{n}\}}.

Remark: The description of the package and its functions are primarily tailored to computing the (complementary) cdf of the two-sided KS statistic, \eqn{D_{n}}.
It should be noted however that one can compute the (complementary) cdf for the one-sided KS statistics \eqn{D_{n}^{-}} or \eqn{D_{n}^{+}} (cf., Dimitrova, Kaishev, Tan (2020)) by appropriately specifying correspondingly \eqn{A_{i} = 0} for all \eqn{i} or \eqn{B_{i} = 1} for all \eqn{i}, in the function \code{\link{ks_c_cdf_Rcpp}}.
}
\details{

The Exact-KS-FFT method underlying \pkg{KSgeneral} is based on expressing the p-value \eqn{P(D_{n} \ge q)} in terms of an appropriate rectangle probability with respect to the uniform order statistics, as noted by Gleser (1985) for \eqn{P(D_{n} > q)}.
The latter representation is used to express \eqn{P(D_{n} \ge q)} via a double-boundary non-crossing probability for a homogeneous Poisson process, with intensity \eqn{n}, which is then efficiently computed using FFT, ensuring total run-time of order \eqn{O(n^{2}log(n))} (see Dimitrova, Kaishev, Tan (2020) and also Moscovich and Nadler (2017) for the special case when \eqn{F(x)} is continuous).

\pkg{KSgeneral} represents an R wrapper of the original C++ code due to Dimitrova, Kaishev, Tan (2020) and based on the C++ code developed by Moscovich and Nadler (2017).
The package includes the functions \code{\link{disc_ks_c_cdf}}, \code{\link{mixed_ks_c_cdf}} and \code{\link{cont_ks_c_cdf}} that compute the complementary cdf \eqn{P(D_n \ge q)}, for a fixed \eqn{q}, \eqn{0 \le q \le 1}, when \eqn{F(x)} is purely discrete, mixed or continuous, respectively.
\pkg{KSgeneral} includes also the functions \code{\link{disc_ks_test}}, \code{\link{mixed_ks_test}} and \code{\link{cont_ks_test}} that compute the p-value \eqn{P(D_{n} \ge d_{n})}, where \eqn{d_{n}} is the value of the KS test statistic computed based on a user provided data sample \eqn{\{x_{1}, ..., x_{n}\}}, when \eqn{F(x)} is purely discrete, mixed or continuous, respectively.


The functions \code{\link{disc_ks_test}} and \code{\link{cont_ks_test}} represent accurate and fast (run time \eqn{O(n^{2}log(n))}) alternatives to the functions \code{\link{ks.test}} from the package \pkg{dgof} and the function \code{\link{ks.test}} from the package \pkg{stat}, which compute p-values of \eqn{P(D_{n} \ge d_{n})}, assuming \eqn{F(x)} is purely discrete or continuous, respectively.

The package also includes the function \code{\link{ks_c_cdf_Rcpp}} which gives the flexibility to compute the complementary cdf (p-value) for the one-sided KS test statistics \eqn{D_{n}^{-}} or \eqn{D_{n}^{+}}.
It also allows for faster computation time and possibly higher accuracy in computing \eqn{P(D_{n} \ge q)}.
%Please report any issue arising or bug in the code to senren.tan@cass.city.ac.uk.

}
\author{
\packageAuthor{KSgeneral}

Maintainer: \packageMaintainer{KSgeneral}
}
\references{
Dimitrina S. Dimitrova, Vladimir K. Kaishev, Senren Tan. (2020) "Computing the Kolmogorov-Smirnov Distribution When the Underlying CDF is Purely Discrete, Mixed or Continuous". Journal of Statistical Software, \bold{95}(10): 1-42. doi:10.18637/jss.v095.i10.

Gleser L.J. (1985). "Exact Power of Goodness-of-Fit Tests of Kolmogorov Type for Discontinuous Distributions". Journal of the American Statistical Association, \bold{80}(392), 954-958.

Moscovich A., Nadler B. (2017). "Fast Calculation of Boundary Crossing Probabilities for Poisson Processes". Statistics and Probability Letters, \bold{123}, 177-182.

}

