\name{dhsic.test}
\alias{dhsic.test}

\title{Independence test based on dHSIC}

\description{
  Hypothesis test for finding statistically significant evidence of
  dependence between several variables. Uses the d-variable Hilbert
  Schmidt independence criterion (dHSIC) as measure of dependence.
  Several types of hypothesis tests are included. The null
  hypothesis (H_0) is that all variables are jointly independent.
}
\usage{
dhsic.test(X, Y, K, alpha = 0.05, method = "permutation",
            kernel = "gaussian", B = 1000, pairwise = FALSE,
            bandwidth = 1, matrix.input = FALSE)
}
\arguments{
  \item{X}{
    either a list of at least two numeric matrices or a single numeric
    matrix. The rows of a matrix correspond to the observations of a
    variable. It is always required that there are an equal number of
    observations for all variables (i.e. all matrices have to have the
    same number of rows). If \code{X} is a single numeric matrix than
    one has to specify the second variable as \code{Y} or set
    \code{matrix.input} to "TRUE". See below for more details.
  }
  \item{Y}{
    a numeric matrix if \code{X} is also a numeric
    matrix and omitted if \code{X} is a list.
  }
  \item{K}{
    a list of the gram matrices corresponding to each variable. If
    \code{K} the following inputs \code{X}, \code{Y}, \code{kernel},
    \code{pairwise}, \code{bandwidth} and \code{matrix.input} will be
    ignored.
  }
  \item{alpha}{
    a numeric value in (0,1) specifying the confidence level of the
    hypothesis test.
  }
  \item{method}{
    a character string specifying the type of hypothesis test used. The
    available options are: "gamma" (gamma approximation based
    test), "permutation" (permutation test (slow)), "bootstrap" (bootstrap
    test (slow)) and "eigenvalue" (eigenvalue based test).
  }
  \item{kernel}{
    a vector of character strings specifying the kernels for each
    variable. There exist two pre-defined kernels: "gaussian" (Gaussian kernel
    with median heuristic as bandwidth) and "discrete" (discrete
    kernel). User defined kernels can also be used by passing the
    function name as a string, which will then be matched using
    \code{\link[base]{match.fun}}. If the length of \code{kernel} is
    smaller than the number of variables the kernel specified in
    \code{kernel[1]} will be used for all variables.
  }
  \item{B}{
    an integer value specifying the number of Monte-Carlo iterations
    made in the permutation and bootstrap test. Only relevant if
    \code{method} is set to "permutation" or to "bootstrap".
  }
  \item{pairwise}{
    a logical value indicating whether one should use HSIC with pairwise
    comparisons instead of dHSIC. Can only be true if there are more than
    two variables.
  }
  \item{bandwidth}{
    a numeric value specifying the size of the bandwidth used for the
    Gaussian kernel. Only used if kernel="gaussian.fixed".
  }
  \item{matrix.input}{
    a boolean. If \code{matrix.input} is "TRUE" the input \code{X} is
    assumed to be a matrix in which the columns correspond to the variables.
  }
}
\details{
  The d-variable Hilbert Schmidt independence criterion is a direct
  extension of the standard Hilbert Schmidt independence criterion
  (HSIC) from two variables to an arbitrary number of
  variables. It is 0 if and only if the variables are jointly independent.
  
  4 different statistical hypothesis tests are implemented all with null hypothesis
  (H_0: \code{X[[1]]},...,\code{X[[d]]} are jointly independent) and alternative hypothesis
  (H_A: \code{X[[1]]},...,\code{X[[d]]} are not jointly independent): 
  1. Permutation test for dHSIC: exact level, slow
  2. Bootstrap test for dHSIC: pointwise asymptotic level and pointwise
  consistent, slow
  3. Gamma approximation based test for dHSIC: only approximate, fast
  4. Eigenvalue based test for dHSIC: pointwise asymptotic level and pointwise
  consistent, medium

  The null hypothesis is rejected if \code{statistic} is strictly
  greater than \code{crit.value}.

  If \code{X} is a list with d matrices, the function tests for joint
  independence of the corresponding d random vectors. If \code{X} is a
  matrix and \code{matrix.input} is "TRUE" the functions tests the
  independence between the columns of \code{X}. If \code{X} is a matrix
  and \code{matrix.input} is "FALSE" then \code{Y} needs to be a matrix,
  too; in this case, the function tests the (pairwise) independence
  between the corresponding two random vectors.
  
  For more details see the references.
}
\value{
A list containing the following components:
\item{statistic}{the value of the test statistic}
\item{crit.value}{critical value of the hypothesis test. The null
  hypothesis (H_0: joint independence) is rejected if \code{statistic}
  is greater than \code{crit.value}.} 
\item{p.value}{p-value of the hypothesis test, i.e. the probability that
  a random version of the test statistic is greater than
  \code{statistic} under the calculated null 
  hypothesis (H_0: joint independence) based on the data.}
\item{time}{numeric vector containing computation times. \code{time[1]} is time
  to compute Gram matrix, \code{time[2]} is time to compute dHSIC and
  \code{time[3]} is the time to compute \code{crit.value} and
  \code{p.value.}}
\item{bandwidth}{bandwidth used during the computation. Only relevant if
Gaussian kernel was used.}
}
\references{
  Gretton, A., K. Fukumizu, C. H. Teo, L. Song, B. Sch{\"o}lkopf and
  A. J. Smola (2007). A kernel statistical test of independence. In
  Advances in Neural Information Processing Systems (pp. 585-592).
  
  Pfister, N., P. B{\"u}hlmann, B. Sch{\"o}lkopf and J. Peters (2017).
  Kernel-based Tests for Joint Independence. To appear in the
  Journal of the Royal Statistical Society, Series B.
}
\author{
Niklas Pfister and Jonas Peters
}


\seealso{
  In order to only compute the test statistic without p-values, use the
  function \code{\link{dhsic}}.
}

\examples{
### pairwise independent but not jointly independent (pairwise HSIC vs dHSIC)
set.seed(0)
x <- matrix(rbinom(100,1,0.5),ncol=1)
y <- matrix(rbinom(100,1,0.5),ncol=1)
z <- matrix(as.numeric((x+y)==1)+rnorm(100),ncol=1)
X <- list(x,y,z)

dhsic.test(X, method="permutation",
           kernel=c("discrete", "discrete", "gaussian"),
           pairwise=TRUE, B=1000)$p.value
dhsic.test(X, method="permutation",
           kernel=c("discrete", "discrete", "gaussian"),
           pairwise=FALSE, B=1000)$p.value

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
\keyword{nonparametric}
