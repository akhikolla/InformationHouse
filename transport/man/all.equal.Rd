\name{all.equal (transport objects)}
\alias{all.equal.pgrid}
\alias{all.equal.pp}
\alias{all.equal.wpp}
%  
\title{
  Methods for Judging Near Equality of Objects of Class pgrid, pp or wpp
}
\description{
  Methods for judging near equality of objects of class pgrid or pp or wpp
}
\usage{
\method{all.equal}{pgrid}(target, current, ...)
\method{all.equal}{pp}(target, current, ...)
\method{all.equal}{wpp}(target, current, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{target, current}{
    the objects of the same class to be compared.
}
  \item{\dots}{
    currently without effect.
}
}
\value{
  Either \code{TRUE} or a vector of \code{\link{mode}} \dQuote{character} describing the differences between target and current.
}
\author{
  Dominic Schuhmacher \email{dschuhm1@uni-goettingen.de}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link[=all.equal]{all.equal (base)}}, \code{\link{compatible}}
}
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
