\name{plot.ivah}

\alias{plot.ivah}

\title{
Plots result of G-estimation in causal AH model. 
}

\description{
This is a \code{plot} method for class \code{"ivah"}. It only supports
objects fitted with \code{estmethod="g"}. 
}

\usage{
\method{plot}{ivah}(x, gof=FALSE, CI.level=0.95, ...)
}

\arguments{
  \item{x}{
an object of class \code{"ivah"}, fitted with \code{estmethod="g"}.
} 
  \item{gof}{
should we plot the goodness-of-fit process? If not, then \eqn{B(t)} is plotted
together with confidence intervals.
}  
  \item{CI.level}{
level for the confidence intervals.
} 
 \item{\dots}{
not used.
}
}

\author{
Arvid Sjolander and Torben Martinussen.
}

\examples{

##See documentation for ivah.

}
