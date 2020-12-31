\name{summary.ivbounds}

\alias{summary.ivbounds}

\title{
Summarizes instrumental variable estimation
}

\description{
This is a \code{summary} method for class \code{"ivbounds"}.
}

\usage{
\method{summary}{ivbounds}(object, ...)
}

\arguments{
  \item{object}{
an object of class \code{"ivbounds"}.
} 
 \item{\dots}{
not used.
}
}

\details{
Provides the lower and and upper bounds for 
\deqn{p_0=p(Y_0=1)}
\deqn{p_1=p(Y_1=1)}
\deqn{\textrm{CRD}=p_1-p_0}
\deqn{\textrm{CRR}=p_1/p_0}
\deqn{\textrm{COR}=\frac{p_1/(1-p_1)}{p_0/(1-p_0)}}

}

\author{
Arvid Sjolander
}

\examples{

##See documentation for ivbounds.

}
