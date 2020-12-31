\name{ks_c_cdf_Rcpp}
\alias{ks_c_cdf_Rcpp}
\title{
R function calling directly the C++ routines that compute the complementary cumulative distribution function of the two-sided (or one-sided, as a special case) Kolmogorov-Smirnov statistic, when the cdf under the null hypothesis is arbitrary (i.e., purely discrete, mixed or continuous)
}
\description{
Function calling directly the C++ routines that compute the complementary cdf for the one-sample two-sided Kolmogorov-Smirnov statistic, given the sample size \code{n} and the file "Boundary_Crossing_Time.txt" in the working directory.
The latter file contains \eqn{A_{i}} and \eqn{B_{i}}, \eqn{i = 1, ..., n}, specified in Steps 1 and 2 of the Exact-KS-FFT method (see Equation (5) in Section 2 of Dimitrova, Kaishev, Tan (2020)).
The latter values form the n-dimensional rectangular region for the uniform order statistics (see Equations (3), (5) and (6) in Dimitrova, Kaishev, Tan (2020)), namely
\eqn{P(D_{n}\ge q) = 1 - P(A_{i} \le U_{(i)} \le B_{i}, 1 \le i \le n) = 1 - P(g(t) \le nU_{n}(t) \le h(t), 0 \le t \le 1)},
where the upper and lower boundary functions \eqn{h(t)}, \eqn{g(t)} are defined as
\eqn{h(t) = \sum_{i=1}^{n}1_{(A_{i} < t)}}, \eqn{g(t) = \sum_{i=1}^{n}1_{(B_{i} \le t)}},
or equivalently, noting that \eqn{h(t)} and \eqn{g(t)} are correspondingly left and right continuous functions, we have
\eqn{\sup\{t\in[0,1]: h(t) < i \} = A_{i}} and \eqn{\inf\{t\in[0,1]: g(t) > i-1 \} = B_{i}}.

Note that on can also compute the (complementary) cdf for the one-sided KS statistics \eqn{D_{n}^{-}} or \eqn{D_{n}^{+}} (cf., Dimitrova, Kaishev, Tan (2020)) by appropriately specifying correspondingly \eqn{A_{i} = 0} for all \eqn{i} or \eqn{B_{i} = 1} for all \eqn{i}, in the function \code{\link{ks_c_cdf_Rcpp}}.
}
\usage{
ks_c_cdf_Rcpp(n)
}
\arguments{
  \item{n}{
  the sample size
}
}
\details{
Note that all calculations here are done directly in C++ and output in R.
That leads to faster computation time, as well as in some cases, possibly higher accuracy (depending on the accuracy of the pre-computed values \eqn{A_{i}} and \eqn{B_{i}}, \eqn{i = 1, ..., n}, provided in the file "Boundary_Crossing_Time.txt") compared to the functions \code{\link{cont_ks_c_cdf}}, \code{\link{disc_ks_c_cdf}}, \code{\link{mixed_ks_c_cdf}}.

Given a random sample \eqn{\{X_{1}, ..., X_{n}\}} of size \code{n} with an empirical cdf \eqn{F_{n}(x)}, the two-sided Kolmogorov-Smirnov goodness-of-fit statistic is defined as \eqn{D_{n} = \sup | F_{n}(x) - F(x) | }, where \eqn{F(x)} is the cdf of a prespecified theoretical distribution under the null hypothesis \eqn{H_{0}}, that \eqn{\{X_{1}, ..., X_{n}\}} comes from \eqn{F(x)}.
The one-sided KS test statistics are correspondingly defined as \eqn{D_{n}^{-} = \sup_{x} (F(x) - F_{n}(x))} and \eqn{D_{n}^{+} = \sup_{x} (F_{n}(x) - F(x))}.

The function \code{\link{ks_c_cdf_Rcpp}} implements the Exact-KS-FFT method, proposed by Dimitrova, Kaishev, Tan (2020), to compute the complementary cdf, \eqn{P(D_{n} \ge q)} at a value \eqn{q}, when \eqn{F(x)} is arbitrary (i.e. purely discrete, mixed or continuous).
It is based on expressing the complementary cdf as
%\eqn{P(D_{n} \ge q) = 1 - P(A_{i} \le U_{(i)} \le B_{i}, 1 \le i \le n)}, where the limites \eqn{a_{i}} and \eqn{b_{i}} in \eqn{F(a_{i}-)} and \eqn{F(b_{i})} are defined as in Step 1 of Dimitrova, Kaishev, Tan (2020).
\eqn{P(D_{n} \ge q) = 1 - P(A_{i} \le U_{(i)} \le B_{i}, 1 \le i \le n)}, where \eqn{A_{i}} and \eqn{B_{i}} are defined as in Step 1 of Dimitrova, Kaishev, Tan (2020).

The complementary cdf is then re-expressed in terms of the conditional probability that a homogeneous Poisson process, \eqn{\xi_{n}(t)} with intensity \eqn{n} will not cross an upper boundary \eqn{h(t)} and a lower boundary \eqn{g(t)}, given that \eqn{\xi_{n}(1) = n} (see Steps 2 and 3 in Section 2.1 of Dimitrova, Kaishev, Tan (2020)). This conditional probability is evaluated using FFT in Step 4 of the method in order to obtain the value of the complementary cdf \eqn{P(D_{n} \ge q)}.
This algorithm ensures a total worst-case run-time of order \eqn{O(n^{2}log(n))} which makes it highly computationally efficient compared to other known algorithms developed for the special cases of continuous or purely discrete \eqn{F(x)}.


The values \eqn{A_{i}} and \eqn{B_{i}}, \eqn{i = 1, ..., n}, specified in Steps 1 and 2 of the Exact-KS-FFT method (see Dimitrova, Kaishev, Tan (2020), Section 2) must be pre-computed (in R or, if needed, using alternative softwares offering high accuracy, e.g. Mathematica) and saved in a file with the name "Boundary_Crossing_Time.txt" (in the current working directory).

The function \code{\link{ks_c_cdf_Rcpp}} is called in R and it first reads the file "Boundary_Crossing_Time.txt" and then computes the value for the complementaty cdf
\eqn{P(D_{n}\ge q) = 1 - P(A_{i} \le U_{(i)} \le B_{i}, 1 \le i \le n) = 1 - P(g(t) \le nU_{n}(t) \le h(t), 0 \le t \le 1)} in C++ and output in R (or as noted above, as a special case, computes the value of the complementary cdf \eqn{P(D_{n}^{+} \ge q) = 1 - P(A_{i} \le U_{(i)} \le 1, 1 \le i \le n)} or \eqn{P(D_{n}^{-} \ge q) = 1 - P(0 \le U_{(i)} \le B_{i}, 1 \le i \le n)}).
}
\value{
Numeric value corresponding to \eqn{P(D_{n}\ge q) = 1 - P(A_{i} \le U_{(i)} \le B_{i}, 1 \le i \le n) = 1 - P(g(t) \le \eta_{n}(t) \le h(t), 0 \le t \le 1)} (or, as a special case, to \eqn{P(D_{n}^{+} \ge q)} or \eqn{P(D_{n}^{-} \ge q)}), given a sample size \code{n} and the file "Boundary_Crossing_Time.txt" containing \eqn{A_{i}} and \eqn{B_{i}}, \eqn{i = 1, ..., n}, specified in Steps 1 and 2 of the Exact-KS-FFT method (see Dimitrova, Kaishev, Tan (2020), Section 2).
}
\references{
Dimitrina S. Dimitrova, Vladimir K. Kaishev, Senren Tan. (2020) "Computing the Kolmogorov-Smirnov Distribution When the Underlying CDF is Purely Discrete, Mixed or Continuous". Journal of Statistical Software, \bold{95}(10): 1-42. doi:10.18637/jss.v095.i10.

Moscovich A., Nadler B. (2017). "Fast Calculation of Boundary Crossing Probabilities for Poisson Processes". Statistics and Probability Letters, \bold{123}, 177-182.

}


\examples{

## Computing the complementary cdf P(D_{n} >= q)
## for n = 10 and q = 0.1, when F(x) is continuous,
## In this case,
## B_i = (i-1)/n + q
## A_i =  i/n - q


n <- 10
q <- 0.1
up_rec <- ((1:n)-1)/n + q
low_rec <- (1:n)/n - q
df <- data.frame(rbind(up_rec, low_rec))
write.table(df,"Boundary_Crossing_Time.txt", sep = ", ",
                  row.names = FALSE, col.names = FALSE)
ks_c_cdf_Rcpp(n)

}
