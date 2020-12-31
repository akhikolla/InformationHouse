\name{dlevene.bi}
\alias{dlevene.bi}
\title{
Calculate Levene's density for a bi-allelic variant
}
\description{
Program \code{dlevene.bi} calculates Levene's density (P(AB|A)) for a
bi-allelic variant.
}
\usage{
dlevene.bi(x)
}
\arguments{
  \item{x}{a vector of genotype counts (AA,AB,BB)}
}
\value{a single real number}
\references{
Levene, H. (1949) On a matching problem arising in genetics. Annals of
Mathematical Statistics, 20,  pp. 91-94.
}
\author{Jan Graffelman (jan.graffelman@upc.edu)}
\seealso{\code{\link{dlevene}},\code{\link{HWExact}}}
\examples{
x <- c(AA=298,AB=489,BB=213)
prob <- dlevene.bi(x)
print(prob)
}
\keyword{distribution}

