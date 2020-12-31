\name{rsp}
\alias{rsp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Distribution-Free Phase I Analysis of Univariate Data based on
  Recursive Segmentation and Permutation
}
\description{
\code{rsp} implements the Phase I method described in Capizzi and Masarotto (2013).
}
\usage{
rsp(y, plot = TRUE, L = 1000, seed = 11642257, alpha = 0.05,
    maxsteps = min(50, round(NROW(y)/15)), lmin = max(5, min(10, round(NROW(y)/10))))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
Phase I data; \code{y} can be either (i) a vector or a 1xm matrix in the case of
individual observations or (ii) a nxm matrix for subgrouped data
(n observations gathered at m time points).
}
  \item{plot}{
    logical; if \code{TRUE}, the diagnostic plot is displayed. 
}
  \item{L}{
    integer; the number of random permutations used to compute the p-values.
}
  \item{seed}{
    positive integer; if not \code{NA}, the RNG's state is resetted
    using \code{seed}. The current \code{.Random.seed} will be preserved.
}
  \item{alpha}{
    real; the significance level used to compute the level and scale
    estimates; if one of the p-values is greater than
    \code{alpha}, the corresponding estimate is a constant.
}
  \item{maxsteps}{
    integer; the maximum number of step shifts which the procedure tries to detect.
}
  \item{lmin}{
    integer; the minimum length of a step.
}
}

\value{
  A list with elements
  \item{p}{the adjusted p-values}
  \item{stat}{the summary statistics (a mx2 matrix)}
  \item{fitted}{the (possibly time-variant) estimates of the process level and
  scale (a mx2 matrix).}
}
\author{
\packageAuthor{dfphase1}.
}
\references{
  G. Capizzi, G. Masarotto (2013),
  \dQuote{Phase I Distribution-Free Analysis of Univariate Data}.
  \emph{Journal of Quality Technology}, \bold{45}, pp. 273-284.
}


\examples{
# Individual observations with a transient level change
set.seed(112233)
level <- c(rep(0,20),rep(3,10),rep(0,20))
x <- level+rt(50,4)
rsp(x)
# Individual observations with a scale step change
scale <- c(rep(1,25),rep(3,25))
x <- scale*rt(50,4)
rsp(x)
data(fe)
rsp(fe)
data(colonscopy)
rsp(colonscopy)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
\keyword{statistical process monitoring}
\keyword{statistical process control}
\keyword{distribution-free methods}
