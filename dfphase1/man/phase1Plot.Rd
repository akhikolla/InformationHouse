\name{phase1Plot}
\alias{phase1Plot}
\alias{mphase1Plot}
\title{
  Plot of Phase 1 data
}
\description{
  \code{phase1Plot} and \code{mphase1Plot}
  plot univariate or multivariate Phase 1 observations,
  organized as required by  the \code{dfphase1} package. 
}
\usage{
phase1Plot(x)

mphase1Plot(x, layout = c(1, p))
}
\arguments{
  \item{x}{
    \code{phase1Plot}: a nxm numeric matrix or a numeric vector of
    length m.

    \code{mphase1Plot}: a pxnxm data numeric array
    or a pxm numeric matrix.

    Here, p denotes the number of variables, n the size
    of each subgroup and m the number of subgroups.
}
  \item{layout}{
    an integer vector describing the multi-panel (and
    possible multi-page) layout. See the third example below.
}
}

\author{
\packageAuthor{dfphase1}.
}
\examples{
  x <- matrix(rt(5*20,5),5)
  x[,10] <- x[,10]+3
  phase1Plot(x)
  # a data set with many variables
  x <- array(rnorm(20*5*50),c(20,5,50))+10*(1:20)
  mphase1Plot(x)
  # it is better to organize the plot on two pages
  if (interactive()) old <- grDevices::devAskNewPage(TRUE)
  mphase1Plot(x,c(2,5,2))
  if (interactive()) grDevices::devAskNewPage(old)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}

