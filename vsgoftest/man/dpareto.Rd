\name{dpareto}
\alias{dpareto}
\alias{ppareto}
\alias{qpareto}
\alias{rpareto}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{The Pareto distribution}
\description{Density, cumulative distribution function, quantile function and random generation for the Pareto distribution.
}
\usage{
dpareto(x, mu, c, log = FALSE)
ppareto(q, mu, c, lower.tail = TRUE, log.p = FALSE)
qpareto(p, mu, c, lower.tail = TRUE, log.p = FALSE)
rpareto(n, mu, c)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x,q}{(\code{numeric, vector}) a vector of quantiles.}
  \item{p}{(\code{numeric, vector}) a vector of probablities.}
  \item{n}{(\code{numeric, vector}) sample size to be generated.}
  \item{mu}{(\code{numeric, single value}) the shape parameter.}
  \item{c}{(\code{numeric, single value}) the scale parameter.}
  \item{log, log.p}{(\code{logical, single value}) if TRUE, probabilities are given as \eqn{\log(p)}. Default is FALSE.}
  \item{lower.tail}{(\code{logical, single value}) if TRUE (default), probabilities are \eqn{P(X\le x)}; otherwise \eqn{P(X>x)}.}
}
\details{The pareto distribution with shape parameter \eqn{\mu>0} and scale parameter \eqn{c>0} has density
\deqn{f(x) = \mu c^\mu x^{-1-\mu}, }
for \eqn{x \ge 0}.
}
\value{
dpareto gives the density, ppareto gives the distribution function, qpareto gives the quantile function, and rpareto generates random deviates.

The length of the result is determined by n for rnorm, and is the maximum of the lengths of the numerical arguments for the other functions.
}
\references{

Arnold, B.C. Pareto distribution, \emph{International Cooperative Publishing House, Fairland} (1983). 

Philbrick, S.W. A practical guide to the single parameter Pareto distribution. \emph{Proceedings of the Casualty Actuarial Society LXXII}, \bold{44}, 44-85 (1985). 

}

\author{J. Lequesne \email{justine.lequesne@unicaen.fr}
}

\examples{
n<- 100
rpareto(n,mu=2,c=1)
 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{distribution}
