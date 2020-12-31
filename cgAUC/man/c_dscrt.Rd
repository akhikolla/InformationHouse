\name{c_dscrt}
\alias{c_dscrt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
c_dscrt}
\description{
discrete function, when variable is discrete.}
\usage{
c_dscrt(y, z, l)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
The potential variables. It is a matrix with column of values of a variables. It should be standardized in this application.}
  \item{z}{
The gold standard variable. It should be standardized.}
  \item{l}{
Linear combination. A vector.}
}
\details{
Discrete function, when variable is discrete.}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{theta.h.p}{The estimate of theta when variable is discrete.}
\item{var}{The variance of estimate of theta.}
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
function(y, z, l) {
    .Call('cgAUC_c_dscrt', PACKAGE = 'cgAUC', y, z, l)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ c_dscrt }
