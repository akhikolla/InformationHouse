\name{estim.mix}
\alias{estim.mix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimation of a mixture model of MCAR and MNAR values in each column of a data matrix.
}
\description{
This function allows estimating a mixture model of MCAR and MNAR values in each column of data sets similar to the ones which can be studied in MS-based quantitative proteomics. Such data matrices contain intensity values of identified peptides.
}
\usage{
estim.mix(tab, tab.imp, conditions, x.step.mod=150, x.step.pi=150,
nb.rei=200)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tab}{A data matrix containing numeric and missing values. Each column of this matrix is assumed to correspond to an experimental sample, and each row to an identified peptide.}
  \item{tab.imp}{A matrix where the missing values of \code{tab} have been imputed under the assumption that they are all MCAR. For instance, such a matrix can be obtained by using the function \code{\link{impute.slsa}} of this package.}
  \item{conditions}{A vector of factors indicating the biological condition to which each column (experimental sample) belongs.}
  \item{x.step.mod}{The number of points in the intervals used for estimating the cumulative distribution functions of the mixing model in each column.}
  \item{x.step.pi}{The number of points in the intervals used for estimating the proportion of MCAR values in each column.}
  \item{nb.rei}{The number of initializations of the minimization algorithm used to estimate the proportion of MCAR values (see Details).}
}
\details{
This function aims to estimate the following mixture model in each column:

\eqn{F_{tot}(x)=\pi_{na}\times F_{na}(x)+(1-\pi_{na})\times F_{obs}(x)}

\eqn{F_{na}(x)=\pi_{mcar}\times F_{tot}(x)+(1-\pi_{mcar})\times F_{mnar}(x)}

where \eqn{\pi_{na}} is the proportion of missing values, \eqn{\pi_{mcar}} is the proportion of MCAR values, \eqn{F_{tot}} is the cumulative distribution function (cdf) of the complete values, \eqn{F_{na}} is the cdf of the missing values, \eqn{F_{obs}} is the cdf of the observed values, and \eqn{F_{mnar}} is the cdf of the MNAR values.

To estimate this model, a first step consists to compute a rough estimate of \eqn{F_{na}}  by assuming that all missing values are MCAR (thanks to the argument \code{tab.imp}). This rough estimate is noted \eqn{\hat{F}_{na}}.

In a second step, the proportion of MCAR values is estimated. To do so, the ratio

\eqn{\hat{\pi}(x)=(1-\hat{F}_{na}(x))/(1-\hat{F}_{tot}(x))}

is computed for different \eqn{x}, where

\eqn{\hat{F}_{tot}(x)=\pi_{na}\times \hat{F}_{na}(x)+(1-\pi_{na})\times \hat{F}_{obs}(x)}

with \eqn{\hat{F}_{obs}} the empirical cdf of the observed values.

Next, the following minimization is performed:

\eqn{\min_{1>k>0,a>0,d>0}f(k,a,d)}

where

\eqn{f(k,a,d)=\sum_x \frac{1}{s(x)^2}\times [\hat{\pi}(x)-k-(1-k)\frac{\exp(-a\times [x-lower_x]^d)}{1-\hat{F}_{tot}(x)}]^2}

where \eqn{s(x)^2} is an estimate of the asymptotic variance of \eqn{\hat{\pi}(x)}, \eqn{lower_x} is an estimate of the minimum of the complete values. To perform this minimization, the function \code{optim} with the method "L-BFGS-B" is used. Because it is function of its initialization, it is possible to reinitialize a number of times the minimisation algorithm with the argument \code{nb.rei}: the parameters leading to the lowest minimum are next kept.

Once \code{k}, \code{a} and \code{d} are estimated, one can use several methods to estimate \eqn{\pi_{mcar}}: it is estimated
by \eqn{k};
}
\value{
A list composed of:
  \item{abs.pi}{A numeric matrix containing the intervals used for estimating the ratio

  \code{(1-F_na(x))/(1-F_tot(x))}

  in each column.}
  \item{pi.init}{A numeric matrix containing the estimated ratios

  \code{(1-F_na(x))/(1-F_tot(x))}

  where \code{x} belongs to \code{abs.pi[,j]} for each sample \code{j}.}
  \item{var.pi.init}{A numeric matrix containing the estimated asymptotic variances of \code{pi.init}.}
  \item{trend.pi.init}{A numeric matrix containing the estimated trend of the model used in the minimization algorithm.}
  \item{abs.mod}{A numeric vector containing the interval used for estimating the mixture models in each column.}
  \item{pi.na}{A numeric vector containing the proportions of missing values in each column.}
  \item{F.na}{A numeric matrix containing the estimated cumulative distribution functions of missing values in each column on the interval \code{abs.mod}.}
  \item{F.tot}{A numeric matrix containing the estimated cumulative distribution functions of complete values in each column on the interval \code{abs.mod}.}
  \item{F.obs}{A numeric matrix containing the estimated cumulative distribution functions of observed values in each column on the interval \code{abs.mod}.}
  \item{pi.mcar}{A numeric vector containing the estimations of the proportion of MCAR values in each column.}
  \item{MinRes}{A numeric matrix containing the three parameters of the model used in the minimization algorithm (three first rows) and the value of minimized function.}
}
\author{
Quentin Giai Gianetto <quentin2g@yahoo.fr>
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{impute.slsa}}
}
\examples{
\donttest{
#Simulating data
res.sim=sim.data(nb.pept=2000,nb.miss=600);

#Imputation of missing values with a MCAR-devoted algorithm: here the slsa algorithm
dat.slsa=impute.slsa(tab=res.sim$dat.obs,conditions=res.sim$condition,repbio=res.sim$repbio);

#Estimation of the mixture model
res=estim.mix(tab=res.sim$dat.obs, tab.imp=dat.slsa, conditions=res.sim$condition);
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Missing value analysis}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
