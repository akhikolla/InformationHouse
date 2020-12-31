\name{robmixglm}
\alias{robmixglm}
\title{Fits a Robust Generalized Linear Model and Variants}
\description{
Fits robust generalized linear models and variants described in Beath (2018).
}
\usage{
robmixglm(formula, family = c("gaussian", "binomial", "poisson", 
"gamma", "truncpoisson", "nbinom"), data, offset = NULL, 
quadpoints = 21, notrials = 50, EMTol = 1.0e-4, cores  =  max(detectCores() - 1,  1), 
verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{formula}{Model formula}
  \item{family}{Distribution of response}
  \item{data}{Data frame from which variables are obtained}
  \item{offset}{Offset to be incorporated in the linear predictor.}
  \item{quadpoints}{Number of quadrature points used in the Gauss-Hermite integration.}
  \item{notrials}{Number of random starting values to be used for EM}
  \item{EMTol}{Relative change in likelihood for completion of EM algorithm before switching to quasi-Newton}
  \item{cores}{Number of cores to be used for parallel evaluation of starting values}
  \item{verbose}{Print out diagnostic information? This includes the likelihood and parameter estimates for each EM run.}
}
\value{
  robmixglm object. This contains
\item{call}{Call to function}
\item{family}{Family of model to be fitted}
  \item{X}{model matrix}
  \item{Y}{response}
  \item{offset}{Offset.}
  \item{mf}{Model frame.}
  \item{mt}{Model terms.}
  \item{xlevels}{Levels for factors.}
  \item{quadpoints}{Number of quadrature points used in the Gauss-Hermite integration.}
  \item{notrials}{Number of random starting values to be used for EM}
  \item{EMTol}{Relative change in likelihood for completion of EM algorithm before switching to quasi-Newton}
\item{verbose}{Was verbose output requested?}
}
\details{
Fits robust generalized models assuming that data is a mixture of standard observations and outlier absorvations, which belong to an overdispersed model. For gaussian, binomial, Poisson, truncated Poisson and gamma, the overdispersed component achieved through including a random effect as part of the linear predictor. For gaussian models this is equivalent to assuming that outlier observations have a higher residual variance, which is specified as $/sigma^2+/tau^2$ to constrain the residual variance to be higher than in the standard model. For negative binomial data the outlier component is also a negative binomial model, but with a lower value of $/theta$ resulting in greater variation.
}
\references{Beath,  K. J. A mixture-based approach to robust analysis of generalised linear models,  Journal of Applied Statistics,  45(12),  2256-2268 (2018) DOI: 10.1080/02664763.2017.1414164}
\author{Ken Beath}
\examples{
\donttest{
library(MASS)
data(forbes)
forbes.robustmix <- robmixglm(100*log10(pres)~bp, data = forbes, cores = 1)
}}
\keyword{multivariate}
