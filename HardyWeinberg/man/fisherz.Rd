\name{fisherz}
\alias{fisherz}
\title{Fisher's z transformation
}
\description{
Calculates Fisher's z transformation for a correlation coefficient
}
\usage{
fisherz(r)
}
\arguments{
  \item{r}{a correlation coefficient
}
}
\value{
 a real number
}
\author{ Jan Graffelman (jan.graffelman@upc.edu) }
\seealso{\code{\link{cor}}}
\examples{
r <- 0.5
print(fisherz(r))
}
\keyword{misc}

