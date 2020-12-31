\name{pi.mcar.logit}
\alias{pi.mcar.logit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimating the proportion of MCAR values in a sample using a logit model.
}
\description{
This function allows estimating the proportion of MCAR values in a sample using a logit model.
}
\usage{
pi.mcar.logit(tab,conditions)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tab}{A data matrix containing numeric and missing values. Each column of this matrix is assumed to correspond to an experimental sample, and each row to an identified peptide.}
  \item{conditions}{A vector of factors indicating the biological condition to which each column (experimental sample) belongs.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A list composed of:
\item{pi.mcar}{The estimated proportion of MCAR values.}
\item{coef1}{The estimated intercept of each logit model estimated in a sample.}
\item{coef2}{The estimated coefficient of each logit model estimated in a sample.}
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Quentin Giai Gianetto <quentin2g@yahoo.fr>
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{estim.mix}}
}
\examples{
\donttest{
#Simulating data
res.sim=sim.data(nb.pept=2000,nb.miss=600);

#Proportion of MCAR values in each sample
pi.mcar.logit(tab=res.sim$dat.obs,conditions=res.sim$conditions)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Missing value analysis}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
