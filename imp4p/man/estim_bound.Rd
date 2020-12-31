\name{estim.bound}
\alias{estim.bound}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimation of lower and upper bounds for missing values.
}
\description{
This function allows estimating lower and upper bounds for missing values of an input matrix. It can be used before to use the functions \code{\link{prob.mcar}} and \code{\link{prob.mcar.tab}}.
}
\usage{
estim.bound(tab, conditions, q=0.95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tab}{A data matrix containing numeric and missing values. Each column of this matrix is assumed to correspond to the intensities measured in an experimental sample, and each row to the ones of an identified peptide.}
  \item{conditions}{A vector of factors indicating the biological condition to which each sample belongs.}
  \item{q}{A numeric value allowing to define confidence intervals for missing values (see Details).}
}
\details{
In each condition, this function estimates lower and upper bounds for missing values of row \code{i} by:

\code{upper(i)=max(tab[i,]);}

\code{lower(i)=max(tab[i,])-quant_diff(q);}

where \code{quant_diff(q)} corresponds to a quantile value of the differences between the maximum and the minimum of the observed values for all the peptides in the condition. As a result, if \code{q} is close to 1, \code{quant_diff(q)} represents an extrem value between the maximum and the minimum of the intensity values in a condition for a peptide.
}
\value{
A list composed of:
  \item{tab.lower}{A matrix with the lower bounds for each missing value in \code{tab}.}
  \item{tab.upper}{A matrix with the upper bounds for each missing value in \code{tab}.}
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
\code{\link{prob.mcar}},\code{\link{prob.mcar.tab}}
}
\examples{
\donttest{
#Simulating data
res.sim=sim.data(nb.pept=2000,nb.miss=600,pi.mcar=0.2,para=0.5,nb.cond=2,nb.repbio=3,
nb.sample=5,m.c=25,sd.c=2,sd.rb=0.5,sd.r=0.2);

data=res.sim$dat.obs;
cond=res.sim$conditions;

#Estimation of lower and upper bounds for each missing value
res=estim.bound(data,conditions=cond);
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Missing value analysis}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
