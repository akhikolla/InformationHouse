%% File Name: immer_hrm.Rd
%% File Version: 0.69

\name{immer_hrm}
\alias{immer_hrm}
\alias{summary.immer_hrm}
\alias{plot.immer_hrm}
\alias{logLik.immer_hrm}
\alias{IRT.likelihood.immer_hrm}
\alias{IRT.posterior.immer_hrm}
\alias{anova.immer_hrm}

%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Hierarchical Rater Model (Patz et al., 2002)
}

\description{
Estimates the hierarchical rater model (HRM; Patz et al., 2002; see Details) with
Markov Chain Monte Carlo using a Metropolis-Hastings algorithm.
}

\usage{
immer_hrm(dat, pid, rater, iter, burnin, N.save=3000, prior=NULL,  est.a=FALSE,
           est.sigma=TRUE,  est.mu=FALSE, est.phi="a", est.psi="a",
           MHprop=NULL, theta_like=seq(-10,10,len=30), sigma_init=1, print_iter=20 )

\method{summary}{immer_hrm}(object, digits=3, file=NULL, ...)

\method{plot}{immer_hrm}(x,...)

\method{logLik}{immer_hrm}(object,...)

\method{anova}{immer_hrm}(object,...)

\method{IRT.likelihood}{immer_hrm}(object,...)

\method{IRT.posterior}{immer_hrm}(object,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
Data frame with item responses
}
  \item{pid}{
Person identifiers
}
  \item{rater}{
Rater identifiers
}
  \item{iter}{
Number of iterations
}
  \item{burnin}{
Number of burnin iterations
}
\item{N.save}{Maximum number of samples to be saved.}
  \item{prior}{
Parameters for prior distributions
}
\item{est.a}{Logical indicating whether \eqn{a} parameter should be estimated.}
\item{est.sigma}{Logical indicating whether \eqn{\sigma} parameter should be
estimated.}
\item{est.mu}{Optional logical indicating whether the mean \eqn{\mu}
of the trait \eqn{\theta} should be estimated.}
\item{est.phi}{Type of \eqn{\phi _{ir}} parameters to be estimated.
    If \code{est.phi="a"}, then \eqn{\phi_{ir}} is estimated for all items and
    all raters. If \code{est.phi="r"}, then \eqn{\phi_{ir}=\phi_r} is
    rater specific, while for \code{est.phi="i"} it is item specific
    (\eqn{\phi_{ir}=\phi_i}). In case of \code{est.phi="n"}, no \eqn{\phi}
    parameters are estimated and all \eqn{\phi} parameters are fixed at 0.
}
\item{est.psi}{Type of \eqn{\psi_{ir}} parameters to be estimated. The
    arguments follow the same conventions as \code{est.phi}, but also
    allows \code{est.psi="e"} (exchangeable) which means \eqn{\psi_{ir}=\psi},
    i.e assuming the same \eqn{\psi} parameter for all items and raters.
}
  \item{MHprop}{
Parameters for Metropolis Hastings sampling. The standard deviation of the
proposal distribution is adaptively computed (Browne & Draper, 2000).
}
\item{theta_like}{Grid of \eqn{\theta} values to be used for likelihood
approximation}
\item{sigma_init}{Initial value for \code{sigma}}
\item{print_iter}{Integer indicating that after each \code{print_iter}th
iteration output on the console should be displayed.}
\item{object}{Object of class \code{immer_hrm}}
\item{digits}{Number of digits after decimal to print}
\item{file}{Name of a file in which the output should be sunk}
\item{x}{Object of class \code{immer_hrm}}
\item{\dots}{Further arguments to be passed. See
\code{\link[sirt:plot.mcmc.sirt]{sirt::plot.mcmc.sirt}} for options in \code{plot}.
}
}

\details{
The hierarchical rater model is defined at the level of persons
\deqn{P( \xi _{pi}=\xi  | \theta_p ) \propto \exp ( \xi
    \cdot a_i \cdot \theta_p - b_{ik}  ) } where \eqn{\theta_p} is normally distributed with mean
\eqn{\mu} and standard deviation \eqn{\sigma}.

At the level of ratings, the model is defined as
\deqn{P( X_{pir}=x | \theta_p, \xi_{pi} ) \propto \exp(
- ( x - \xi_{pi} - \phi_{ir} )^2 / ( 2 \cdot \psi_{ir} ) ) }
}

\value{
A list with following entries
  \item{person}{Data frame with estimated person parameters}
  \item{item}{Data frame with estimated item parameters}
  \item{rater_pars}{Data frame with estimated rater parameters}
  \item{est_pars}{Estimated item and trait distribution parameters arranged
         in vectors and matrices.}
  \item{sigma}{Estimated standard deviation \eqn{\sigma} of trait \eqn{\theta}}
  \item{mu}{Estimated mean \eqn{\mu} of trait \eqn{\theta}}
  \item{mcmcobj}{Object of class \code{mcmc.list} for \pkg{coda} package.}
  \item{summary.mcmcobj}{Summary of all parameters}
  \item{EAP.rel}{EAP reliability}
  \item{ic}{Parameters for information criteria}
  \item{f.yi.qk}{Individual likelihood evaluated at \code{theta_like}}
  \item{f.qk.yi}{Individual posterior evaluated at \code{theta_like}}
  \item{theta_like}{Grid of \eqn{\theta} values for likelihood approximation}
  \item{pi.k}{Discretized \eqn{\theta} distribution}
   \item{like}{Log-likelihood value}
  \item{MHprop}{Updated parameters in Metropolis-Hastings sampling}
}


\references{
Browne, W. J., & Draper, D. (2000). Implementation and performance issues in
the Bayesian and likelihood fitting of multilevel models.
\emph{Computational Statistics, 15}, 391-420.

Patz, R. J., Junker, B. W., Johnson, M. S., & Mariano, L. T. (2002).
The hierarchical rater model for rated test items and its
application to large-scale educational assessment data.
\emph{Journal of Educational and Behavioral Statistics, 27}(4), 341-384.
}

%% Browne, W. J., & Draper, D. (2006). A comparison of Bayesian and
%% likelihood-based methods for fitting multilevel models.
%% \emph{Bayesian Analysis}, \bold{1(3)}, 473-514.



%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
Simulate the HRM with \code{\link{immer_hrm_simulate}}.
}


\examples{
\dontrun{
library(sirt)
library(TAM)

#############################################################################
# EXAMPLE 1: Simulated data using the immer::immer_hrm_simulate() function
#############################################################################

# define data generating parameters
set.seed(1997)
N <- 500  # number of persons
I <- 4    # number of items
R <- 3    # number of raters
K <- 3    # maximum score
sigma <- 2  # standard deviation
theta <- stats::rnorm( N, sd=sigma )  # abilities
# item intercepts
b <- outer( seq( -1.5, 1.5, len=I), seq( -2, 2, len=K), "+" )
# item loadings
a <- rep(1,I)
# rater severity parameters
phi <- matrix( c(-.3, -.2, .5), nrow=I, ncol=R, byrow=TRUE )
phi <- phi + stats::rnorm( phi, sd=.3 )
phi <- phi - rowMeans(phi)
# rater variability parameters
psi <- matrix( c(.1, .4, .8), nrow=I, ncol=R, byrow=TRUE )
# simulate HRM data
data <- immer::immer_hrm_simulate( theta, a, b, phi=phi, psi=psi )
pid <- data$pid
rater <- data$rater
dat <- data[, - c(1:2) ]

#----------------------------------------------------------------
#*** Model 1: estimate HRM with equal item slopes
iter <- 3000
burnin <- 500
mod1 <- immer::immer_hrm( dat, pid, rater, iter=iter,  burnin=burnin )
summary(mod1)
plot(mod1,layout=2,ask=TRUE)

# relations among convergence diagnostic statistics
par(mfrow=c(1,2))
plot( mod1$summary.mcmcobj$PercVarRatio, log(mod1$summary.mcmcobj$effSize), pch=16)
plot( mod1$summary.mcmcobj$PercVarRatio, mod1$summary.mcmcobj$Rhat, pch=16)
par(mfrow=c(1,1))

# extract individual likelihood
lmod1 <- IRT.likelihood(mod1)
str(lmod1)
# extract log-likelihood value
logLik(mod1)

# write coda files into working directory
sirt::mcmclist2coda(mod1$mcmcobj, name="hrm_mod1")

#----------------------------------------------------------------
#*** Model 2: estimate HRM with estimated item slopes
mod2 <- immer::immer_hrm( dat, pid, rater, iter=iter,  burnin=burnin,
            est.a=TRUE, est.sigma=FALSE)
summary(mod2)
plot(mod2,layout=2,ask=TRUE)

# model comparison
anova( mod1, mod2 )

#----------------------------------------------------------------
#*** Model 3: Some prior specifications
prior <- list()
# prior on mu
prior$mu$M <- .7
prior$mu$SD <- 5
# fixed item parameters for first item
prior$b$M <-  matrix( 0, nrow=4, ncol=3 )
prior$b$M[1,] <- c(-2,0,2)
prior$b$SD <-  matrix( 10, nrow=4, ncol=3 )
prior$b$SD[1,] <- 1E-4
# estimate model
mod3 <- immer::immer_hrm( dat, pid, rater, iter=iter,  burnin=burnin, prior=prior)
summary(mod3)
plot(mod3)

#----------------------------------------------------------------
#*** Model 4: Multi-faceted Rasch models in TAM package

# create facets object
facets <- data.frame( "rater"=rater )

#-- Model 4a: only main rater effects
form <- ~ item*step + rater
mod4a <- TAM::tam.mml.mfr( dat, pid=pid, facets=facets, formulaA=form)
summary(mod4a)

#-- Model 4b: item specific rater effects
form <- ~ item*step + item*rater
mod4b <- TAM::tam.mml.mfr( dat, pid=pid, facets=facets, formulaA=form)
summary(mod4b)

#----------------------------------------------------------------
#*** Model 5: Faceted Rasch models with sirt::rm.facets

#--- Model 5a: Faceted Rasch model with only main rater effects
mod5a <- sirt::rm.facets( dat, pid=pid, rater=rater )
summary(mod5a)

#--- Model 5b: allow rater slopes for different rater discriminations
mod5b <- sirt::rm.facets( dat, pid=pid, rater=rater, est.a.rater=TRUE )
summary(mod5b)

#############################################################################
# EXAMPLE 2: data.ratings1 (sirt package)
#############################################################################

data(data.ratings1, package="sirt")
dat <- data.ratings1

# set number of iterations and burnin iterations
set.seed(87)
iter <- 1000
burnin <- 500
# estimate model
mod <- immer::immer_hrm( dat[,  paste0("k",1:5) ], pid=dat$idstud, rater=dat$rater,
             iter=iter, burnin=burnin )
summary(mod)
plot(mod, layout=1, ask=TRUE)
plot(mod, layout=2, ask=TRUE)
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Hierarchical rater model}
\keyword{plot}
\keyword{summary}
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
