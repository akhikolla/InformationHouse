\name{extractAIC}
\title{Extract AIC from a Fitted Model}
%
\alias{extractAIC.robmixglm}
\usage{
\method{extractAIC}{robmixglm}(fit,  scale,  k  =  2,  \dots)
}
\arguments{
  \item{fit}{fitted \code{robmixglm} model.}
 \item{scale}{ignored.}
 \item{k}{numeric specifying the \sQuote{weight} of the
   \emph{equivalent degrees of freedom} (\eqn{\equiv}{ = :} \code{edf})
   part in the AIC formula.}
 \item{\dots}{further arguments (currently unused).}
}
\description{
  Computes the (generalized) AIC for a fitted \code{robmixglm} model. Used in \code{step},  otherwise use \code{AIC}.
}

\author{Ken Beath}
\seealso{
  \code{\link{extractAIC}}, \code{\link{step}}
}
\examples{\donttest{
library(MASS)
data(forbes)
forbes.robustmix <- robmixglm(bp~pres, data = MASS::forbes, cores = 1)
extractAIC(forbes.robustmix)
}}
\keyword{models}
