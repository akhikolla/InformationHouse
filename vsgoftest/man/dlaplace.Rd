\name{dlaplace}
\alias{dlaplace}
\alias{plaplace}
\alias{qlaplace}
\alias{rlaplace}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{The Laplace distribution}
\description{Density, cumulative distribution function, quantile function and random generation for the laplace distribution.
}
\usage{
dlaplace(x, mu, b, log = FALSE)
plaplace(q, mu, b, lower.tail = TRUE, log.p = FALSE)
qlaplace(p, mu, b, lower.tail = TRUE, log.p = FALSE)
rlaplace(n, mu, b)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x,q}{(\code{numeric, vector}) a vector of quantiles.}
  \item{p}{(\code{numeric, vector}) a vector of probablities.}
  \item{n}{(\code{numeric, vector}) sample size to be generated.}
  \item{mu}{(\code{numeric, single value}) the location parameter.}
  \item{b}{(\code{numeric, single value}) the scale parameter.}
  \item{log, log.p}{(\code{logical, single value}) if TRUE, probabilities are given as \eqn{\log(p)}. Default is FALSE.}
  \item{lower.tail}{(\code{logical, single value}) if TRUE (default), probabilities are \eqn{P(X\le x)}; otherwise \eqn{P(X>x)}.}
}
\details{The laplace distribution with shape parameter \eqn{\mu>0} and scale parameter \eqn{b>0} has density
\deqn{f(x) = \frac{1}{2b}\exp(-\frac{|x-\mu|}{b}),  \quad x \in R.}
}
\value{
dlaplace gives the density, plaplace gives the distribution function, qlaplace gives the quantile function, and rlaplace generates random deviates.

The length of the result is determined by n for rnorm, and is the maximum of the lengths of the numerical arguments for the other functions.
}

\author{J. Lequesne \email{justine.lequesne@unicaen.fr}
}

\examples{
set.seed(1)
rlaplace(100,mu=2,b=1)
}
\keyword{distribution}
