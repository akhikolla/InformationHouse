%% File Name: immer_ccml.Rd
%% File Version: 0.16

\name{immer_ccml}
\alias{immer_ccml}
\alias{summary.immer_ccml}
\alias{coef.immer_ccml}
\alias{vcov.immer_ccml}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Composite Conditional Maximum Likelihood Estimation for the Partial Credit Model
with a Design Matrix for Item Parameters
}


\description{
Estimates the partial credit model with a design matrix for item
parameters with composite conditional maximum likelihood estimation.
The estimation uses pairs of items \eqn{X_i} and \eqn{X_j} and considers conditional
likelihoods \eqn{P(X_i=k, X_j=h | \theta) / P( X_i + X_j=k+h| \theta )}. By using
this strategy, the trait \eqn{\theta} cancels out (like in conditional maximum
likelihood estimation). The proposed strategy is a generalization of the Zwinderman (1995)
composite conditional maximum likelihood approach of the Rasch model to the
partial credit model. See Varin, Reid and Firth (2011) for a general introduction to
composite conditional maximum likelihood estimation.
}


\usage{
immer_ccml( dat, weights=NULL, irtmodel="PCM", A=NULL, b_fixed=NULL, control=NULL )

\method{summary}{immer_ccml}(object, digits=3, file=NULL, ...)

\method{coef}{immer_ccml}(object, ...)

\method{vcov}{immer_ccml}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
Data frame with polytomous item responses \eqn{0,1,\ldots, K}
}
  \item{weights}{
Optional vector of sampling weights
}
\item{irtmodel}{Model string for specifying the item response model}
\item{A}{
Design matrix (items \eqn{\times} categories \eqn{\times} basis parameters).
Entries for categories are for \eqn{1,\ldots,K}
}
  \item{b_fixed}{
Matrix with fixed \eqn{b} parameters
}
\item{control}{Control arguments for optimization function
\code{\link[stats:nlminb]{stats::nlminb}}
}
\item{object}{Object of class \code{immer_ccml}}
\item{digits}{Number of digits after decimal to print}
\item{file}{Name of a file in which the output should be sunk}
\item{\dots}{Further arguments to be passed.}
}




\details{
The function estimates the partial credit model as
\eqn{P(X_i=h | \theta ) \propto \exp( h \theta - b_{ih} )} with
\eqn{b_{ih}=\sum_l a_{ihl} \xi_l} where the values \eqn{a_{ihl}}
are included in the design matrix \code{A} and \eqn{\xi_l} denotes
basis item parameters.
}

\value{
List with following entries (selection)

\item{coef}{Item parameters}
\item{vcov}{Covariance matrix for item parameters}
\item{se}{Standard errors for item parameters}
\item{nlminb_result}{Output from optimization with
\code{\link[stats:nlminb]{stats::nlminb}}}
\item{suff_stat}{Used sufficient statistics}
\item{ic}{Information criteria}

}

\references{
Varin, C., Reid, N., & Firth, D. (2011). An overview of composite likelihood methods.
\emph{Statistica Sinica, 21}, 5-42.

Zwinderman, A. H. (1995). Pairwise parameter estimation in Rasch models.
\emph{Applied Psychological Measurement, 19}(4), 369-375.
}

%\author{
%%  ~~who you are~~
%}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See \code{\link[sirt:rasch.pairwise.itemcluster]{sirt::rasch.pairwise.itemcluster}}
of an implementation of the composite conditional maximum likelihood approach for the
Rasch model.
}


\examples{
#############################################################################
# EXAMPLE 1: Partial credit model with CCML estimation
#############################################################################

library(TAM)

data(data.gpcm, package="TAM")
dat <- data.gpcm

#-- initial MML estimation in TAM to create a design matrix
mod1a <- TAM::tam.mml(dat, irtmodel="PCM2")
summary(mod1a)

#* define design matrix
A <- - mod1a$A[,-1,-1]
A <- A[,,-1]
str(A)

#-- estimate model
mod1b <- immer::immer_ccml( dat, A=A)
summary(mod1b)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
