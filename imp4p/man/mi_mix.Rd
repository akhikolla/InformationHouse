\name{mi.mix}
\alias{mi.mix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multiple imputation from a matrix of probabilities of being MCAR for each missing value.
}
\description{
This function allows imputing data sets with a multiple imputation strategy.
}
\usage{
mi.mix(tab, tab.imp, prob.MCAR, conditions, repbio=NULL, reptech=NULL, nb.iter=3, nknn=15,
weight=1, selec="all", siz=500, ind.comp=1, methodMCAR="mle", q=0.95,
progress.bar=TRUE, details=FALSE, ncp.max=5, maxiter = 10, ntree = 100,
variablewise = FALSE, decreasing = FALSE, verbose = FALSE, mtry = floor(sqrt(ncol(tab))),
replace = TRUE,classwt = NULL, cutoff = NULL, strata = NULL, sampsize = NULL,
nodesize = NULL, maxnodes = NULL,xtrue = NA, parallelize = c('no', 'variables',
'forests'), methodMNAR="igcda",q.min = 0.025, q.norm = 3, eps = 0, distribution = "unif",
param1 = 3, param2 = 1, R.q.min=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tab}{A data matrix containing numeric and missing values. Each column of this matrix is assumed to correspond to an experimental sample, and each row to an identified peptide.}
  \item{tab.imp}{A matrix where the missing values of \code{tab} have been imputed under the assumption that they are all MCAR. For instance, such a matrix can be obtained from the function \code{\link{impute.slsa}} of this package.}
  \item{prob.MCAR}{A matrix of probabilities that each missing value is MCAR. For instance such a matrix can be obtained from the function \code{\link{prob.mcar.tab}} of this package.}
  \item{conditions}{A vector of factors indicating the biological condition to which each column (experimental sample) belongs.}
  \item{repbio}{A vector of factors indicating the biological replicate to which each column belongs. Default is NULL (no experimental design is considered).}
  \item{reptech}{A vector of factors indicating the technical replicate to which each column belongs. Default is NULL (no experimental design is considered).}
  \item{nb.iter}{The number of iterations used for the multiple imputation method.}
  \item{nknn}{The number of nearest neighbours used in the SLSA algorithm (see \code{\link{impute.slsa}}).}
  \item{selec}{A parameter to select a part of the dataset to find nearest neighbours between rows. This can be useful for big data sets (see \code{\link{impute.slsa}}).}
  \item{siz}{A parameter to select a part of the dataset to perform imputations with a MCAR-devoted algorithm. This can be useful for big data sets. Note that \code{siz} needs to be inferior to \code{selec}.}
  \item{weight}{The way of weighting in the algorithm (see \code{\link{impute.slsa}}).}
  \item{ind.comp}{If \code{ind.comp=1}, only nearest neighbours without missing values are selected to fit linear models (see \code{\link{impute.slsa}}). Else, they can contain missing values.}
  \item{methodMCAR}{The method used for imputing MCAR data. If \code{methodi="mle"} (default), then the \code{\link{impute.mle}} function is used (imputation using an EM algorithm). If \code{methodi="pca"}, then the \code{\link{impute.PCA}} function is used (imputation using Principal Component Analysis). If \code{methodi="rf"}, then the \code{\link{impute.RF}} function is used (imputation using Random Forest). Else, the \code{\link{impute.slsa}} function is used (imputation using Least Squares on nearest neighbours).}
  \item{methodMNAR}{The method used for imputing MNAR data. If \code{methodMNAR="igcda"} (default), then the \code{\link{impute.igcda}} function is used. Else, the \code{\link{impute.pa}} function is used.}
  \item{q}{A quantile value (see \code{\link{impute.igcda}}).}
  \item{progress.bar}{If \code{TRUE}, a progress bar is displayed.}
  \item{details}{If \code{TRUE}, the function gives a list of three values: \code{imputed.matrix} a matrix with the average of imputed values for each missing value, \code{sd.imputed.matrix} a matrix with the standard deviations of imputed values for each missing value, \code{all.imputed.matrices} an array with all the \code{nb.iter} matrices of imputed values that have been generated.}
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
At each iteration, a matrix indicating the MCAR values is generated by Bernouilli distributions having parameters given by the matrix \code{prob.MCAR}. The generated MCAR values are next imputed thanks to the matrix \code{tab.imp}. For each row containing MNAR values, the other rows are imputed thanks to the function  \code{\link{impute.igcda}} and, next, the considered row is imputed thanks to one of the MCAR-devoted imputation methods (\code{\link{impute.mle}}, \code{\link{impute.RF}}, \code{\link{impute.PCA}} or \code{\link{impute.slsa}}). So, the function \code{\link{impute.igcda}} allows to deform the correlation structure of the dataset in view to be closer to that of the true values, while the MCAR-devoted imputation method will impute by taking into account this modified correlation structure.
}
\value{
The input matrix \code{tab} with average imputed values instead of missing values if \code{details=FALSE} (default). If \code{details=TRUE}, a list of three values: \code{imputed.matrix} a matrix with the average of imputed values for each missing value, \code{sd.imputed.matrix} a matrix with the standard deviations of imputed values for each missing value, \code{all.imputed.matrices} an array with all the \code{nb.iter} matrices of imputed values that have been generated.
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
res.sim=sim.data(nb.pept=5000,nb.miss=1000);

#Fast imputation of missing values with the impute.rand algorithm
dat.rand=impute.rand(tab=res.sim$dat.obs,conditions=res.sim$condition);

#Estimation of the mixture model
res=estim.mix(tab=res.sim$dat.obs, tab.imp=dat.rand, conditions=res.sim$condition);

#Computing probabilities to be MCAR
born=estim.bound(tab=res.sim$dat.obs,conditions=res.sim$condition);
proba=prob.mcar.tab(tab.u=born$tab.upper,res=res);

#Multiple imputation strategy with 3 iterations (can be time consuming in function of the data set!)
data.mi=mi.mix(tab=res.sim$dat.obs, tab.imp=dat.rand, prob.MCAR=proba, conditions=
res.sim$conditions, repbio=res.sim$repbio, nb.iter=3);
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Imputation methods}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
