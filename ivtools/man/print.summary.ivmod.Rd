\name{print.summary.ivmod}

\alias{print.summary.ivmod}

\title{
Prints summary of instrumental variable estimation
}

\description{
This is a \code{print} method for class \code{"summary.ivmod"}.
}

\usage{
\method{print}{summary.ivmod}(x, digits=max(3L, getOption("digits")-3L), 
  signif.stars=getOption("show.signif.stars"), ...)
}

\arguments{
  \item{x}{
an object of class \code{"summary.ivmod"}.
}
  \item{digits}{
the number of significant digits to use.
}
  \item{signif.stars}{
logical. If TRUE, "significance stars" are printed for each coefficient.
}
  \item{\dots}{
not used.
}
}

\author{
Arvid Sjolander
}


\examples{

##See documentation for ivglm, ivcoxph and ivah.

}
