\name{impute.igcda}
\alias{impute.igcda}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Imputing missing values by assuming that the distribution of complete values is Gaussian in each column of an input matrix. This algorithm is named "Imputation under a Gaussian Complete Data Assumption" (IGCDA).
}
\description{
This function allows imputing missing values under the assumption that the distribution of complete values has to be Gaussian in each column.

Note that the imputed values are not necessary small values (compared to observed values).
}
\usage{
impute.igcda(tab, tab.imp, conditions, q=0.95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tab}{A numeric vector or matrix with observed and missing values.}
  \item{tab.imp}{A matrix where the missing values of \code{tab} have been imputed under the assumption that they are all MCAR. For instance, such a matrix can be obtained by using the function \code{\link{impute.slsa}} of this package.}
  \item{conditions}{A vector of factors indicating the biological condition to which each column (experimental sample) belongs.}
  \item{q}{A quantile value (see Details).}
}
\details{
The mean and variance of the Gaussian distribution are determined using a linear regression between the quantiles of the observed values \code{q_{obs}} and the ones of the standard normal distribution \code{q_{N(0,1)}}.

The quantile value is used for determining the minimum of imputed values. This minimum is determined by the minimum observed value in the dataset minus \code{quant_diff(q)} where \code{quant_diff(q)} corresponds to a quantile value of the differences between the maximum and the minimum of the observed values for all the peptides in the condition. As a result, if \code{q} is close to 1, \code{quant_diff(q)} represents an extrem value between the maximum and the minimum of the intensity values in a condition for a peptide.
}
\value{
The numeric input matrix with imputed values. The distribution of the intensity values in each of its columns is supposed to be Gaussian.
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

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
\donttest{
#Simulating data
res.sim=sim.data(nb.pept=2000,nb.miss=600);

#Imputation of missing values with a MCAR-devoted algorithm: here the slsa algorithm
dat.slsa=impute.slsa(tab=res.sim$dat.obs,conditions=res.sim$condition,repbio=res.sim$repbio);

#Imputation of missing values under a Gaussian assumption
dat.gauss=impute.igcda(tab=res.sim$dat.obs, tab.imp=dat.slsa, conditions=res.sim$conditions);
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Imputation methods}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
