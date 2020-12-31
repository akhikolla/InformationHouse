\name{ifisherz}
\alias{ifisherz}
\title{Inverse Fisher z transformation
}
\description{
Calculates the inverse of Fisher's z transformation
}
\usage{
ifisherz(y)
}
\arguments{
  \item{y}{a real number}
}
\value{
  a correlation coefficient in the range (-1,1)
}
\author{ Jan Graffelman (jan.graffelman@upc.edu) }
\seealso{\code{\link{cor}}}
\examples{
r <- 0.5
print(ifisherz(fisherz(r)))
}
\keyword{misc}

