%% File Name: immer_latent_regression.Rd
%% File Version: 0.12

\name{immer_latent_regression}
\alias{immer_latent_regression}
\alias{summary.immer_latent_regression}
\alias{coef.immer_latent_regression}
\alias{vcov.immer_latent_regression}
\alias{logLik.immer_latent_regression}
\alias{anova.immer_latent_regression}


%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Unidimensional Latent Regression
}

\description{
Fits a unidimensional latent regression \eqn{\theta_{ig}=Y_{ig} \bold{\beta} + \varepsilon_{ig}}
with group-specific variances \eqn{Var(\varepsilon _{ig} )=\sigma^2_g} based on
the individual likelihood of a fitted model.
}

\usage{
immer_latent_regression(like, theta=NULL, Y=NULL, group=NULL, weights=NULL,
   conv=1e-05, maxit=200, verbose=TRUE)

\method{summary}{immer_latent_regression}(object, digits=3, file=NULL, ...)

\method{coef}{immer_latent_regression}(object, ...)

\method{vcov}{immer_latent_regression}(object, ...)

\method{logLik}{immer_latent_regression}(object, ...)

\method{anova}{immer_latent_regression}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{like}{
Matrix containing the individual likelihood \eqn{L( \bold{X} | \theta )}
}
  \item{theta}{
Grid of \eqn{\bold} values
}
  \item{Y}{
Predictor matrix
}
  \item{group}{
Group identifiers
}
  \item{weights}{
Optional vector of weights
}
  \item{conv}{
Convergence criterion
}
  \item{maxit}{
Maximum number of iterations
}
  \item{verbose}{
Logical indicating whether progress should be displayed
}
\item{object}{Object of class \code{immer_latent_regression}}
\item{digits}{Number of digits after decimal to print}
\item{file}{Name of a file in which the output should be sunk}
\item{\dots}{Further arguments to be passed.}
}

%%\details{
%% describe some details
%%}

\value{
List containing values (selection)

\item{coef}{Parameter vector}
\item{vcov}{Covariance matrix for estimated parameters}
\item{beta}{Regression coefficients}
\item{gamma}{Standard deviations}
\item{beta_stat}{Data frame with \eqn{\bold{\beta}} parameters}
\item{gamma_stat}{Data frame with standard deviations}
\item{ic}{Information criteria}
\item{deviance}{Deviance}
\item{N}{Number of persons}
\item{G}{Number of groups}
\item{group}{Group identifier}
\item{iter}{Number of iterations}

}

\references{
Adams, R. J., & Wu, M. L. (2007). The mixed-coefficients multinomial logit model.
A generalized form of the Rasch model. In M. von Davier & C. H. Carstensen (Eds.):
\emph{Multivariate and mixture distribution Rasch models: Extensions and applications} (pp. 55-76).
New York: Springer.
}


%\author{
%%  ~~who you are~~
%}

\note{
The \code{IRT.likelihood} method can be used for
extracting the individual likelihood.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See \code{\link[TAM:tam.latreg]{TAM::tam.latreg}} for latent regression estimation
in the \pkg{TAM} package.
}


\examples{
\dontrun{
#############################################################################
# EXAMPLE 1: Latent regression for Rasch model with simulated data
#############################################################################

library(sirt)

#-- simulate data
set.seed(9877)
I <- 15  # number of items
N <- 700 # number of persons per group
G <- 3   # number of groups
b <- seq(-2,2,len=I)
group <- rep( 1:G, each=N)
mu <- seq(0,1, length=G)
sigma <- seq(1, 1.5, length=G)
dat <- sirt::sim.raschtype( stats::rnorm( N*G, mean=mu[group], sd=sigma[group] ), b)

#-- estimate Rasch model with JML
mod1 <- immer::immer_jml( dat )
summary(mod1)

#-- compute individual likelihood
like1 <- IRT.likelihood(mod1)

#-- estimate latent regression
mod2 <- immer::immer_latent_regression( like=like1, group=group)
summary(mod2)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
