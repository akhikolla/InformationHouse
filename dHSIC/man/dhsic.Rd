\name{dhsic}
\alias{dhsic}

\title{d-variable Hilbert Schmidt independence criterion - dHSIC}

\description{
  The d-variable Hilbert Schmidt independence criterion (dHSIC) is
  a non-parametric measure of dependence between an arbitrary number of
  variables. In the large sample limit the value of dHSIC is 0 if the
  variables are jointly independent and positive if there is a
  dependence. It is therefore able to detect any type of dependence
  given a sufficient amount of data.
}
\usage{
dhsic(X, Y, K, kernel = "gaussian", bandwidth = 1, matrix.input = FALSE)
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
    \code{K} specified the other inputs will have no effect on the
    computations.
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
  variables. It is 0 if and only if all the variables are jointly
  independent. This function computes an estimator of dHSIC, which
  converges to the actual dHSIC in the large sample limit. It is
  therefore possible to detect any type of dependence in the large
  sample limit.

  If \code{X} is a list with d matrices, the function computes dHSIC for
  the corresponding d random vectors. If \code{X} is a
  matrix and \code{matrix.input} is "TRUE" the functions dHSIC for the
  columns of \code{X}. If \code{X} is a matrix and \code{matrix.input}
  is "FALSE" then \code{Y} needs to be a matrix, too; in this case, the
  function computes the dHSIC (HSIC) for the corresponding two random vectors.
  
  For more details see the references.
}
\value{
A list containing the following components:
\item{dHSIC}{the value of the empirical estimator of dHSIC}
\item{time}{numeric vector containing computation times. \code{time[1]}
  is time to compute Gram matrix and \code{time[2]} is time to compute
  dHSIC.}
\item{bandwidth}{bandwidth used during computations. Only relevant if
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
  In order to perform hypothesis tests based on dHSIC use the function \code{\link{dhsic.test}}.
}

\examples{
### Three different input methods
set.seed(0)
x <- matrix(rnorm(200),ncol=2)
y <- matrix(rbinom(100,30,0.1),ncol=1)
# compute dHSIC of x and y (x is taken as a single variable)
dhsic(list(x,y),kernel=c("gaussian","discrete"))$dHSIC
dhsic(x,y,kernel=c("gaussian","discrete"))$dHSIC
# compute dHSIC of x[,1], x[,2] and y
dhsic(cbind(x,y),kernel=c("gaussian","discrete"), matrix.input=TRUE)$dHSIC

### Using a user-defined kernel (here: sigmoid kernel)
set.seed(0)
x <- matrix(rnorm(500),ncol=1)
y <- x^2+0.02*matrix(rnorm(500),ncol=1)
sigmoid <- function(x_1,x_2){
  return(tanh(sum(x_1*x_2)))
}
dhsic(x,y,kernel="sigmoid")$dHSIC

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{nonparametric}
