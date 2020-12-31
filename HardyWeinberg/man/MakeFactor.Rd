\name{MakeFactor}
\alias{MakeFactor}
\title{Make factors from genotyping data}
\description{
  \code{MakeFactor} converts bi-allelic genetic marker data, whether coded numerically
  as (0,1,2) or as (GG,GT,TT), etc. into standard factors coded as AA,
  AB, BB.
}
\usage{
MakeFactor(x, coding = c(0, 1, 2))
}
\arguments{
  \item{x}{A vector containing genotyping results}
  \item{coding}{Describes the numerical coding of the genotype data in
    order AA, AB and BB. Only relevant if \code{x} is numerical}
}
\details{If \code{x} is a factor, it will be coerced to a factor with
  levels AA, AB and BB. Important detail: the produced factors will have
  only those levels that are observed in the data. E.g., if genotyping
  results only consist of (0,1), then the resulting factor will not have
  level BB (which would be an empty category)
}
\value{A factor variable}
\author{ Jan Graffelman \email{jan.graffelman@upc.edu} }
\seealso{
\code{\link{MakeCounts}}
}
\examples{
y <- c(1,1,0,0,2,2)
data.frame(y,MakeFactor(y))

y <- c(2,2,3,3,1,1)
data.frame(y,MakeFactor(y,coding=c(1,2,3)))

data(Markers)
data.frame(Markers[,1],MakeFactor(Markers[,1],coding=c(1,2,3)))
}
\keyword{misc}

