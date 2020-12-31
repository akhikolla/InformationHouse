%% File Name: immer_jml.Rd
%% File Version: 0.49

\name{immer_jml}
\alias{immer_jml}
\alias{summary.immer_jml}
\alias{IRT.likelihood.immer_jml}
\alias{logLik.immer_jml}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Joint Maximum Likelihood Estimation for the Partial Credit Model
with a Design Matrix for Item Parameters
and \eqn{\varepsilon}-Adjustment Bias Correction
}


\description{
Estimates the partial credit model with a design matrix for item
parameters with joint maximum likelihood (JML). The  \eqn{\varepsilon}-adjustment
bias correction is implemented with reduces bias of the
JML estimation method (Bertoli-Barsotti, Lando & Punzo, 2014).
}

\usage{
immer_jml(dat, A=NULL, maxK=NULL, center_theta=TRUE, b_fixed=NULL, irtmodel="PCM",
     pid=NULL, rater=NULL, eps=0.3, est_method="eps_adj", maxiter=1000,
     conv=1e-05, max_incr=3, maxiter_update=10, maxiter_line_search=6,
     conv_update=1e-05, verbose=TRUE, use_Rcpp=TRUE, shortcut=TRUE)

\method{summary}{immer_jml}(object, digits=3, file=NULL, ...)

\method{logLik}{immer_jml}(object, ...)

\method{IRT.likelihood}{immer_jml}(object, theta=seq(-9,9,len=41), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
Data frame with polytomous item responses \eqn{0,1,\ldots, K}
}
  \item{A}{
Design matrix (items \eqn{\times} categories \eqn{\times} basis parameters).
Entries for categories are for \eqn{1,\ldots,K}
}
\item{maxK}{Optional vector with maximum category per item}
  \item{center_theta}{
Logical indicating whether the trait estimates should be centered
}
  \item{b_fixed}{
Matrix with fixed \eqn{b} parameters
}
\item{irtmodel}{Specified item response model. Can be one of the two
partial credit model parametrizations \code{PCM} and \code{PCM2}.
}
\item{pid}{Person identifier}
\item{rater}{Optional rater identifier}
  \item{eps}{
Adjustment parameter \eqn{\varepsilon}
}
  \item{est_method}{
Estimation method. Can be \code{'eps_adj'} for the \eqn{\varepsilon}-adjustment,
\code{'jml'} for the JML without bias correction and
\code{'jml_bc'} for JML with bias correction.
}
  \item{maxiter}{
Maximum number of iterations
}
  \item{conv}{
Convergence criterion
}
  \item{max_incr}{
Maximum increment
}
  \item{maxiter_update}{
Maximum number of iterations for parameter updates
}
\item{maxiter_line_search}{Maximum number of iterations within line search}
  \item{conv_update}{
Convergence criterion for updates
}
  \item{verbose}{
Logical indicating whether convergence progress should be displayed
}
  \item{use_Rcpp}{
Logical indicating whether \pkg{Rcpp} package should be used
for computation.
}
\item{shortcut}{Logical indicating whether a computational shortcut should be used
for efficiency reasons}
\item{object}{Object of class \code{immer_jml}}
\item{digits}{Number of digits after decimal to print}
\item{file}{Name of a file in which the output should be sunk}
\item{theta}{Grid of \eqn{\theta} values}
\item{\dots}{Further arguments to be passed.}
}




\details{
The function uses the partial credit model as
\eqn{P(X_i=h | \theta ) \propto \exp( h \theta - b_{ih} )} with
\eqn{b_{ih}=\sum_l a_{ihl} \xi_l} where the values \eqn{a_{ihl}}
are included in the design matrix \code{A} and \eqn{\xi_l} denotes
basis item parameters.

The adjustment parameter \eqn{\varepsilon} is applied to the sum score
as the sufficient statistic for the person parameter. In more detail,
extreme scores \eqn{S_p=0} (minimum score) or \eqn{S_p=M_p} (maximum score)
are adjusted to \eqn{S_p^\ast=\varepsilon} or \eqn{S_p^\ast=M_p - \varepsilon},
respectively. Therefore, the adjustment possesses more influence on
parameter estimation for datasets with a small number of items.
}

\value{
List with following entries

\item{b}{Item parameters \eqn{b_{ih}}}
\item{theta}{Person parameters}
\item{theta_se}{Standard errors for person parameters}
\item{xsi}{Basis parameters}
\item{xsi_se}{Standard errors for bias parameters}
\item{probs}{Predicted item response probabilities}
\item{person}{Data frame with person scores}
\item{dat_score}{Scoring matrix}
\item{score_pers}{Sufficient statistics for persons}
\item{score_items}{Sufficient statistics for items}
\item{loglike}{Log-likelihood value}

}

\references{
Bertoli-Barsotti, L., Lando, T., & Punzo, A. (2014). Estimating a Rasch Model via
fuzzy empirical probability functions. In D. Vicari, A. Okada, G. Ragozini &
C. Weihs (Eds.).
\emph{Analysis and Modeling of Complex Data in Behavioral and Social Sciences},
Springer.
}

%\author{
%%  ~~who you are~~
%}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See \code{\link[TAM:tam.jml]{TAM::tam.jml}} for
joint maximum likelihood estimation. The \eqn{varepsilon}-adjustment
is also implemented in \code{\link[sirt:mle.pcm.group]{sirt::mle.pcm.group}}.
}


\examples{
#############################################################################
# EXAMPLE 1: Rasch model
#############################################################################

data(data.read, package="sirt")
dat <- data.read

#---  Model 1: Rasch model with JML and epsilon-adjustment
mod1a <- immer::immer_jml(dat)
summary(mod1a)

\dontrun{
#- JML estimation, only handling extreme scores
mod1b <- immer::immer_jml( dat, est_method="jml")
summary(mod1b)

#- JML estimation with (I-1)/I bias correction
mod1c <- immer::immer_jml( dat, est_method="jml_bc" )
summary(mod1c)

# compare different estimators
round( cbind( eps=mod1a$xsi, JML=mod1b$xsi, BC=mod1c$xsi ), 2 )

#---  Model 2: LLTM by defining a design matrix for item difficulties
A <- array(0, dim=c(12,1,3) )
A[1:4,1,1] <- 1
A[5:8,1,2] <- 1
A[9:12,1,3] <- 1

mod2 <- immer::immer_jml(dat, A=A)
summary(mod2)

#############################################################################
# EXAMPLE 2: Partial credit model
#############################################################################

library(TAM)
data(data.gpcm, package="TAM")
dat <- data.gpcm

#-- JML estimation in TAM
mod0 <- TAM::tam.jml(resp=dat, bias=FALSE)
summary(mod0)

# extract design matrix
A <- mod0$A
A <- A[,-1,]

#-- JML estimation
mod1 <- immer::immer_jml(dat, A=A, est_method="jml")
summary(mod1)

#-- JML estimation with epsilon-adjusted bias correction
mod2 <- immer::immer_jml(dat, A=A, est_method="eps_adj")
summary(mod2)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
