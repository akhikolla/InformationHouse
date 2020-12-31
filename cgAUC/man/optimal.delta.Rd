\name{optimal.delta}
\alias{optimal.delta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
optimal.delta}
\description{
Find the optimal delta.
}
\usage{
optimal.delta(y, z, l, h, ind.d.l)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
The potential variables. It is a matrix with column of values of a variables. It should be standardized in this application.}
  \item{z}{
The gold standard variable. It should be standardized.}
  \item{l}{
Linear combination. A vector.}
  \item{h}{
The value of h falls into (n^(-1/2), n^(-1/5)).}
  \item{ind.d.l}{
Void}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{delta.star}{Optimal delta.}
}
\references{
Chang, YCI. Maximizing an ROC type measure via linear combination of markers when the gold reference is continuous. Statistics in Medicine 2012.\cr
Obuchowski NA. An ROC-type measure of diagnostic accuracy when the gold standard is continuous-scale. Statistics in Medicine 2006; 25:481--493.\cr
Obuchowski N. Estimating and comparing diagnostic tests accuracy when the gold standard is not binary. Statistics in Medicine 2005; 20:3261--3278.\cr
Friedman JH, Popescu BE. Gradient directed regularization for linear regression and classification. Technical Report, Department of Statistics, Stanford University, 2004.\cr
}
\author{
Yu-chia Chang}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (y, z, l, h, ind.d.l) 
{
    l.i = matrix(rep(l, times = 50), nrow = 50, byrow = TRUE)
    delta = seq(0, 5, length = 50)
    m = delta \%*\% t(ind.d.l)
    l.i = l.i + m
    l.i.max = apply(l.i, 1, max)
    l.i = l.i/l.i.max
    theta = rep(0, 50)
    for (i in 2:50) {
        theta[i] = cntin(y, z, l.i[i, ], h)$theta.sh.h.p
    }
    delta.star = delta[which(theta == max(theta))]
    return(delta.star)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ optimal.delta }
