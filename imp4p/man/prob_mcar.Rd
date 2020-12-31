\name{prob.mcar}
\alias{prob.mcar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimation of a vector of probabilities that missing values are MCAR.
}
\description{
This function returns a vector of probabilities that each missing value is MCAR from specified confidence intervals.
}
\usage{
prob.mcar(b.u,absc,pi.na,pi.mcar,F.tot,F.obs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{b.u}{A numeric vector of upper bounds for missing values.}
  \item{absc}{The interval on which is estimated the MCAR data mechanism.}
  \item{pi.na}{The estimated proportion of missing values.}
  \item{pi.mcar}{The estimated proportion of MCAR values among missing values.}
  \item{F.tot}{An estimation of the cumulative distribution function of the complete values on the interval \code{absc}.}
  \item{F.obs}{An estimation of the cumulative distribution function of the missing values on the interval \code{absc}.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A numeric vector of estimated probabilities to be MCAR for missing values assuming upper bounds for them (\code{b.u}). The input arguments \code{absc}, \code{pi.mcar}, \code{pi.na}, \code{F.tot} and \code{F.obs} can be estimated thanks to the function \code{\link{estim.mix}}.
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
#Simulating data
res.sim=sim.data(nb.pept=2000,nb.miss=600);

#Imputation of missing values with the slsa algorithm
dat.slsa=impute.slsa(tab=res.sim$dat.obs,conditions=res.sim$condition,repbio=res.sim$repbio);

#Estimation of the mixture model
res=estim.mix(tab=res.sim$dat.obs, tab.imp=dat.slsa, conditions=res.sim$condition);

#Computing probabilities to be MCAR
born=estim.bound(tab=res.sim$dat.obs,conditions=res.sim$condition);

#Computing probabilities to be MCAR in the first column of result$tab.mod
proba=prob.mcar(b.u=born$tab.upper[,1],absc=res$abs.mod,pi.na=res$pi.na[1],
pi.mcar=res$pi.mcar[1], F.tot=res$F.tot[,1], F.obs=res$F.obs[,1]);
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Missing value analysis}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
