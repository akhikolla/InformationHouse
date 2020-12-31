\name{compatible}
\alias{compatible}
\alias{compatible.pgrid}
\alias{compatible.pp}
\alias{compatible.wpp}
%
\title{
  Test whether Two Objects are Compatible
}
\description{
  Test whether two objects of the same class are \sQuote{of similar shape} so that
  the function \code{\link{transport}} can be applied.
}
\usage{
compatible(target, current, ...)
\method{compatible}{pgrid}(target, current, ...)
\method{compatible}{pp}(target, current, ...)
\method{compatible}{wpp}(target, current, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{target, current}{
    to objects of the same class to be compared.
}
  \item{\dots}{
    currently without effect.
}
}
\value{
  Logical.
}
\author{
  Dominic Schuhmacher \email{dschuhm1@uni-goettingen.de}
}

\seealso{
  \code{\link[=all.equal.pgrid]{all.equal}}
}
