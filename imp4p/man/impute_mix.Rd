\name{impute.mix}
\alias{impute.mix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Imputation using a decision rule under an assumption of a mixture of MCAR and MNAR values.
}
\description{
This function allows imputing data sets with a MCAR-devoted algorithm and a MNAR-devoted algorithm using probabilities that missing values are MCAR. If such a probability is superior to a chosen threshold, then the MCAR-devoted algorithm is used, otherwise it is the MNAR-devoted algorithm.
}
\usage{
impute.mix(tab, prob.MCAR, threshold, conditions, repbio=NULL, reptech=NULL,
methodMCAR="mle",nknn=15,weight=1, selec="all", ind.comp=1, progress.bar=TRUE, q=0.95,
ncp.max=5, maxiter = 10, ntree = 100, variablewise = FALSE, decreasing = FALSE,
verbose = FALSE, mtry = floor(sqrt(ncol(tab))), replace = TRUE,classwt = NULL,
cutoff = NULL, strata = NULL, sampsize = NULL, nodesize = NULL, maxnodes = NULL,
xtrue = NA, parallelize = c('no', 'variables', 'forests'),
methodMNAR="igcda", q.min = 0.025, q.norm = 3, eps = 0, distribution = "unif",
param1 = 3, param2 = 1, R.q.min=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tab}{A data matrix containing numeric and missing values. Each column of this matrix is assumed to correspond to an experimental sample, and each row to an identified peptide.}
  \item{prob.MCAR}{A matrix of probabilities that each missing value is MCAR. For instance such a matrix can be obtained from the function \code{\link{prob.mcar.tab}} of this package.}
  \item{threshold}{A value such that if the probability that a missing value is MCAR is superior to it, then a MCAR-devoted algorithm is used, otherwise it is a MNAR-devoted algorithm that is used.}
  \item{conditions}{A vector of factors indicating the biological condition to which each column (experimental sample) belongs.}
  \item{repbio}{A vector of factors indicating the biological replicate to which each column belongs. Default is NULL (no experimental design is considered).}
  \item{reptech}{A vector of factors indicating the technical replicate to which each column belongs. Default is NULL (no experimental design is considered).}
  \item{methodMCAR}{The method used for imputing MCAR data. If \code{methodi="mle"} (default), then the \code{\link{impute.mle}} function is used (imputation using an EM algorithm). If \code{methodi="pca"}, then the \code{\link{impute.PCA}} function is used (imputation using Principal Component Analysis). If \code{methodi="rf"}, then the \code{\link{impute.RF}} function is used (imputation using Random Forest). Else, the \code{\link{impute.slsa}} function is used (imputation using Least Squares on nearest neighbours).}
    \item{methodMNAR}{The method used for imputing MNAR data. If \code{methodMNAR="igcda"} (default), then the \code{\link{impute.igcda}} function is used. Else, the \code{\link{impute.pa}} function is used.}
  \item{nknn}{The number of nearest neighbours used in the SLSA algorithm (see \code{\link{impute.slsa}}).}
  \item{weight}{The way of weighting in the algorithm (see \code{\link{impute.slsa}}).}
  \item{selec}{A parameter to select a part of the dataset to find nearest neighbours between rows. This can be useful for big data sets (see \code{\link{impute.slsa}}).}
  \item{ind.comp}{If \code{ind.comp=1}, only nearest neighbours without missing values are selected to fit linear models (see \code{\link{impute.slsa}}). Else, they can contain missing values.}
  \item{progress.bar}{If \code{TRUE}, a progress bar is displayed.}
  \item{q}{A quantile value (see \code{\link{impute.igcda}}).}
  \item{ncp.max}{parameter of the \code{\link{impute.PCA}} function.}
  \item{maxiter}{parameter of the \code{\link{impute.RF}} function.}
  \item{ntree}{parameter of the \code{\link{impute.RF}} function.}
  \item{variablewise}{parameter of the \code{\link{impute.RF}} function.}
  \item{decreasing}{parameter of the \code{\link{impute.RF}} function.}
  \item{verbose}{parameter of the \code{\link{impute.RF}} function.}
  \item{mtry}{parameter of the \code{\link{impute.RF}} function.}
  \item{replace}{parameter of the \code{\link{impute.RF}} function.}
  \item{classwt}{parameter of the \code{\link{impute.RF}} function.}
  \item{cutoff}{parameter of the \code{\link{impute.RF}} function.}
  \item{strata}{parameter of the \code{\link{impute.RF}} function.}
  \item{sampsize}{parameter of the \code{\link{impute.RF}} function.}
  \item{nodesize}{parameter of the \code{\link{impute.RF}} function.}
  \item{maxnodes}{parameter of the \code{\link{impute.RF}} function.}
  \item{xtrue}{parameter of the \code{\link{impute.RF}} function.}
  \item{parallelize}{parameter of the \code{\link{impute.RF}} function.}
    \item{q.min}{parameter of the \code{\link{impute.pa}} function.}
  \item{q.norm}{parameter of the \code{\link{impute.pa}} function.}
  \item{eps}{parameter of the \code{\link{impute.pa}} function.}
  \item{distribution}{parameter of the \code{\link{impute.pa}} function.}
  \item{param1}{parameter of the \code{\link{impute.pa}} function.}
  \item{param2}{parameter of the \code{\link{impute.pa}} function.}
  \item{R.q.min}{parameter of the \code{\link{impute.pa}} function.}
}
\details{
The missing values for which \code{prob.MCAR} is superior to a chosen threshold are imputed with one of the MCAR-devoted imputation methods (\code{\link{impute.mle}}, \code{\link{impute.RF}}, \code{\link{impute.PCA}} or \code{\link{impute.slsa}}). The other missing values are considered MNAR and imputed with \code{\link{impute.igcda}}. More details and explanations can be bound in Giai Gianetto (2020).
}
\value{
The input matrix \code{tab} with imputed values instead of missing values.
}
\references{
Giai Gianetto, Q., Wieczorek S., Couté Y., Burger, T. (2020).  A peptide-level multiple imputation strategy accounting for the different natures of missing values in proteomics data. bioRxiv 2020.05.29.122770; doi: \href{https://doi.org/10.1101/2020.05.29.122770}{10.1101/2020.05.29.122770}
}
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

#Fast imputation of missing values with the impute.rand algorithm
dat.rand=impute.rand(tab=res.sim$dat.obs,conditions=res.sim$condition);

#Estimation of the mixture model
res=estim.mix(tab=res.sim$dat.obs, tab.imp=dat.rand, conditions=res.sim$condition);

#Computing probabilities to be MCAR
born=estim.bound(tab=res.sim$dat.obs,conditions=res.sim$condition);
proba=prob.mcar.tab(born$tab.upper,res);

#Imputation under the assumption of MCAR and MNAR values
tabi=impute.mix(tab=res.sim$dat.obs, prob.MCAR=proba, threshold=0.5, conditions=res.sim$conditions,
repbio=res.sim$repbio, methodMCAR="slsa", methodMNAR="igcda", nknn=15, weight=1, selec="all",
ind.comp=1, progress.bar=TRUE);
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Imputation methods}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
