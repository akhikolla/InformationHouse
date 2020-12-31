\name{GenerateSamples}
\alias{GenerateSamples}
\title{Generate genotypic compositions}
\description{
  \code{GenerateSamples} generates all possible genotypic compositions
  (AA,AB,BB) for a given sample size \code{n}.
}
\usage{
GenerateSamples(n = 5)
}
\arguments{
  \item{n}{the desired sample size}
}
\value{returns a matrix with in each row a possible genotypic compostion
  for the given sample size.
}
\author{ Jan Graffelman \email{jan.graffelman@upc.edu} }
\examples{
GenerateSamples(5)
}
\keyword{datagen}

