\name{miss.total.process}
\alias{miss.total.process}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimating the missing data mechanism in a sample.
}
\description{
This function allows estimating the missing data mechanism, i.e. the probability to be missing  in function of the intensity level, from an estimation of a mixture model of MNAR and MCAR values (see \code{\link{estim.mix}} function).
}
\usage{
miss.total.process(abs,pi_na,F_na,F_tot)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{abs}{The interval on which is estimated the missing data mechanism.}
    \item{pi_na}{The proportion of missing values.}
  \item{F_na}{An estimation of the cumulative distribution function of the missing values on the interval \code{abs}.}
\item{F_tot}{An estimation of the cumulative distribution function of the complete values on the interval \code{abs}.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A list composed of:
\item{abs}{The interval on which is estimated the missing data mechanism.}
\item{p}{The estimated probability to be missing in function of the intensity level.}
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

#Imputation of missing values with the slsa algorithm
dat.slsa=impute.slsa(tab=res.sim$dat.obs,conditions=res.sim$condition,repbio=res.sim$repbio);

#Estimation of the mixture model
res=estim.mix(tab=res.sim$dat.obs, tab.imp=dat.slsa, conditions=res.sim$condition);

#Estimating the missing mechanism in the first replicate
mtp=miss.total.process(res$abs.mod,res$pi.na[1],res$F.na[,1],res$F.tot[,1])
plot(mtp$abs,mtp$p,ty="l",xlab="Intensity values",ylab="Estimated probability to be missing")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Missing value analysis}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
