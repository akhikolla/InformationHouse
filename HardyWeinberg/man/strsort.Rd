\name{strsort}
\alias{strsort}
\title{Sort tokens of a set of strings}
\description{
Function \code{strsort} collapses all tokens of a vector of strings in a
single string with sorted tokens}
\usage{
strsort(s)
}
\arguments{
  \item{s}{a vector of character strings}
}
\value{a string}
\author{ Jan Graffelman \email{jan.graffelman@upc.edu} }
\seealso{\code{\link{alleles}}}
\examples{
x <- c("AA","AB","BB","AC","CC")
print(strsort(x))
}
\keyword{misc}