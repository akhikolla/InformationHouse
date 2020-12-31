\name{c_s_h}
\alias{c_s_h}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
c_s_h}
\description{
Smooth function.}
\usage{
c_s_h(t, h)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t}{
A value, the difference between any two subjects.}
  \item{h}{
The value of h falls into (n^(-1/2), n^(-1/5)).}
}
\details{
Smooth function.}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{s_h}{The value of smooth function.}
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
function(t, h) {
    .Call('cgAUC_c_s_h', PACKAGE = 'cgAUC', t, h)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ c_s_h }
