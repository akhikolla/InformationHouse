% $Id: norWood2004.Rd 82 2010-10-24 20:00:48Z Lars $
\name{norWood2004}
\alias{norWood2004}
\docType{data}

\title{
Data: Forestry in Norway
}

\description{
A data set for 113 farmers in forestry in Norway.
}

\usage{data(norWood2004)}

\format{
  A data frame with 113 observations on the following 7 variables.
  \describe{
    \item{\code{firm}}{firm number}
    \item{\code{m}}{Variable cost}
    \item{\code{x}}{Woodland, value of forrest and land}
    \item{\code{y}}{Profit}
    \item{\code{z1}}{Secondary income from ordinary farming}
    \item{\code{z3}}{Age of forrest owner}
    \item{\code{z6}}{Whether there is a long-term plan =1 or not =0}
  }
}


\details{
Collected from farmers in forestry.
}


\note{ Data as \code{.csv} are loaded by the command \code{data} using
\code{ read.table(..., header=TRUE, sep=";")} such that this file
is a semicolon separated file and not a comma separated file.  }


\source{
 Norwegian Agricultural Economics Research Institute.
}


\examples{
data(norWood2004)
## maybe str(norWood2004) ; plot(norWood2004) ...
}

\keyword{datasets}
