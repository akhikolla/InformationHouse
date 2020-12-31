\name{bs.design}
\alias{bs.design}
\concept{bs.design}
\title{Generate a B-spline design matrix}
\description{The function generates a B-spline design matrix with equidistant knots for given degree of the splines and number of basis functions.
}

\usage{

bs.design(x, xl, xr, spline.degree, nbasis, comp = NULL)} 
    
\arguments{
  \item{x}{the positions where spline to be evaluated.}
  \item{xl}{lower intervall boundary where spline functions are relevant.}
  \item{xr}{upper intervall boundary where spline functions are relevant.}
  \item{spline.degree}{(polynomial) degree of the B-splines.}
  \item{nbasis}{number of basis functions used.}
  \item{comp}{Specify if only specific columns of the B-spline design matrix should be returned. Default is NULL and the whole B-spline design matrix is returned.}
}

\value{
  The B-spline design matrix is returned.
}

\author{
Andreas Groll \email{groll@math.lmu.de}
}

\seealso{
  \code{\link{pencoxfrail}}
}

\examples{
x <- rnorm(100)
B <- bs.design(x=x, xl=min(x), xr=max(x), spline.degree=3, nbasis=5)
}
