\name{c_cntin}
\alias{c_cntin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
c_cntin}
\description{
Continue function, when variable was continue.}
\usage{
c_cntin(y, z, l, h)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{The potential variables. It is a matrix with column of values of a variables. It should be standardized in this application.
}
  \item{z}{
The gold standard variable. It should be standardized.}
  \item{l}{
Linear combination. A vector.}
  \item{h}{
The value of h falls into (n^(-1/2), n^(-1/5)).}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{theta.sh.h.p}{The estimate of the theta of Chang(2012).}
\item{var}{The variance of estimate of the theta of Chang(2012).}

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
function(y, z, l, h) {
    .Call('cgAUC_c_cntin', PACKAGE = 'cgAUC', y, z, l, h)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ c_cntin }
