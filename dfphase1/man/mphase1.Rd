\name{mphase1}
\alias{mphase1}
\title{
Distribution-free Phase I analysis of multivariate data 
}
\description{
  Retrospective change point detection using the method  described by
  Capizzi and Masarotto (2016).
}
\usage{
mphase1(x, plot = TRUE, post.signal = TRUE, isolated = dim(x)[2] > 1, step = TRUE,
        alpha = 0.05, gamma = 0.5, K = min(50, round(sqrt(dim(x)[3]))),
        lmin = 5, L = 1000, seed = 11642257)
}
\arguments{
  \item{x}{
    a pxnxm array containing the observations; \code{x[r,j,i]} is
    the jth observation on the rth variable of the ith subgroup. 
}
  \item{plot}{
    logical; if \code{FALSE} the diagnostic plot is not displayed.
}
  \item{post.signal}{
    logical; if \code{FALSE} the diagnostic LASSO-based analysis
    is not performed.
  }
  \item{isolated}{
    logical; if \code{FALSE} isolated shifts are not detected.
  }
  \item{step}{
    logical; if \code{FALSE} step shifts are not detected.
  }
  \item{alpha}{
    real; the acceptable false alarm probability; if the
    observed p-value is greater than \code{alpha}, then
    the estimated mean function is a constant.
}
  \item{gamma}{
    real; the extra penalization for the extended BIC criteria.
}
  \item{K}{
    integer; the maximum number of shifts which the procedure
    tries to detect.
}
  \item{lmin}{
    integer; the minimum length of a step shift.
}
  \item{L}{
    integer; the number of random permutations used to compute the p-values.
}
\item{seed}{
  integer; if not \code{NA}, the RNG's state is re-setted
  using \code{seed}. The current \code{.Random.seed} will be preserved.
}
}
\value{
  Functions \code{mphase1} returns an object of class \code{mphase1} containing
  \item{p.value}{The p-value.}
  \item{Wobs}{The overall test statistic.}
  \item{alasso}{A data-frame containing the result of the
    post-signal diagnosis analysis,i.e., the times and types of shifts
    and the involved variables identified using the adaptive LASSO.}
  \item{forward}{A data frame containing the result of the forward search analysis, i.e.,
  the times and types of the possible shifts as well as
  the elementary test statistics and the estimates of their
  (conditional) means and standard deviations.}
  \item{center, scatter}{The location vector and dispersion matrix used
    to standardize the original data.}
  \item{signed.ranks}{A pxnxm array containing the signed ranks.}
  \item{fitted, residuals}{Two pxnxm arrays containing the fitted
    means and the residuals, i.e., the difference between the
    observations and the fitted values.
  }
}
\references{
    G. Capizzi and G. Masarotto (2016), Phase I Distribution-Free
  Analysis of Multivariate Data, \emph{Technometrics},
  to appear, \url{http://dx.doi.org/10.1080/00401706.2016.1272494}.
}
\author{
\packageAuthor{dfphase1}.
}
\seealso{
 \code{\link{postsignal}}.
}
\examples{
  # A simulated in-control data from a Student's t distribution
  # with 3 degrees of freedom
  set.seed(123)
  x <- sweep(array(rnorm(5*5*50),c(5,5,50)),c(2,3),sqrt(rchisq(5*50,3)/3),"/")
  mphase1(x)
  # Reproduction of the two examples given in Capizzi and Masarotto (2016)
  data(ryan)
  mphase1(ryan)
  data(gravel)
  mphase1(gravel)
}
\keyword{htest}
\keyword{multivariate}
\keyword{statistical process monitoring}
\keyword{statistical process control}
\keyword{distribution-free methods}

