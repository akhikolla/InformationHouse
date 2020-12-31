\name{randfun}
\alias{randfun}
\title{
       Rand Index
}
\description{
Measures to compare the dissimilarity of exact cluster labels (memberships) and estimated cluster labels (memberships) of the observations.
}
\usage{
randfun(lvel, dv)
}
\arguments{
  \item{lvel}{exact cluster labels of the observations }
  \item{dv}{estimated cluster labels of the observations}
}
\value{
a single value between 0 and 1
}
\references{
William M Rand (1971). Objective criteria for the evaluation of clustering methods, \emph{Journal of the American Statistical association}, 66(336):846-850, doi:10.1080/01621459.1971.10482356.
}
\author{ Biplab Paul, Shyamal K. De and Anil K. Ghosh

Maintainer: Biplab Paul<biplab.paul@niser.ac.in>
}
\examples{
   # Measures of dissimilarity:
   exl <- c(rep(0,5), rep(1,5), rep(2,5), rep(3,5))
   el <- c(0,0,1,0,0,1,2,1,0,1,2,2,3,2,2,3,2,3,1,3)
   randfun(exl,el)

   ## outputs:
   #[1] 0.2368421
}
\keyword{ randfun function}