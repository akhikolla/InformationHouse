%% File Name: immer_opcat.Rd
%% File Version: 0.20

\name{immer_opcat}
\alias{immer_opcat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimation of Integer Item Discriminations
}
\description{
Estimates integer item discriminations like in the one-parameter
logistic model (OPLM; Verhelst & Glas, 1995).
See Verhelst, Verstralen and Eggen (1991) for computational
details.
}
\usage{
immer_opcat(a, hmean, min=1, max=10, maxiter=200)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a}{
Vector of estimated item discriminations
}
  \item{hmean}{
Prespecified harmonic mean
}
  \item{min}{
Minimum integer item discrimination
}
  \item{max}{
Maximum integer item discrimination
}
  \item{maxiter}{
Maximum number of iterations
}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{
Vector containing integer item discriminations
}

\references{
Verhelst, N. D. &, Glas, C. A. W. (1995). The one-parameter logistic model.
In G. H. Fischer & I. W. Molenaar (Eds.). \emph{Rasch Models} (pp. 215--238).
New York: Springer.

Verhelst, N. D., Verstralen, H. H. F. M., & Eggen, T. H. J. M. (1991).
\emph{Finding starting values for the item parameters and suitable discrimination
indices in the one-parameter logistic model}.
CITO Measurement and Research Department Reports, 91-10.
}




%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See \code{\link{immer_cml}} for using \code{immer_opcat} to estimate
the one-parameter logistic model.
}


\examples{
#############################################################################
# EXAMPLE 1: Estimating integer item discriminations for dichotomous data
#############################################################################

library(sirt)
data(data.read, package="sirt")
dat <- data.read
I <- ncol(dat)

#--- estimate 2PL model
mod <- sirt::rasch.mml2( dat, est.a=1:I, mmliter=30)
summary(mod)
a <- mod$item$a        # extract (non-integer) item discriminations

#--- estimate integer item discriminations under different conditions
a1 <- immer::immer_opcat( a, hmean=3, min=1, max=6 )
table(a1)
a2 <- immer::immer_opcat( a, hmean=2, min=1, max=3 )
a3 <- immer::immer_opcat( a, hmean=1.5, min=1, max=2 )
#--- compare results
cbind( a, a1, a2, a3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Integer item discriminations}
\keyword{OPLM}
