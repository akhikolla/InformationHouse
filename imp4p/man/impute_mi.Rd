\name{impute.mi}
\alias{impute.mi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Imputation of data sets containing peptide intensities with a multiple imputation strategy.
}
\description{
This function allows imputing data sets containing peptide intensities with a multiple imputation strategy.
}
\usage{
impute.mi(tab, conditions, repbio=NULL, reptech=NULL, nb.iter=3, nknn=15, selec=1000,
siz=900, weight=1, ind.comp=1, progress.bar=TRUE, x.step.mod=300, x.step.pi=300,
nb.rei=100, q=0.95, methodMCAR="mle",ncp.max=5,
maxiter = 10, ntree = 100, variablewise = FALSE,
decreasing = FALSE, verbose = FALSE, mtry = floor(sqrt(ncol(tab))),
replace = TRUE, classwt = NULL, cutoff = NULL, strata = NULL, sampsize = NULL,
nodesize = NULL, maxnodes = NULL, xtrue = NA,
parallelize = c('no', 'variables', 'forests'),
methodMNAR="igcda",q.min = 0.025, q.norm = 3, eps = 0,
distribution = "unif", param1 = 3, param2 = 1, R.q.min=1);
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tab}{A data matrix containing only numeric and missing values. Each column of this matrix is assumed to correspond to an experimental sample, and each row to an identified peptide.}
  \item{conditions}{A vector of factors indicating the biological condition to which each column (experimental sample) belongs.}
  \item{repbio}{A vector of factors indicating the biological replicate to which each column belongs. Default is NULL (no experimental design is considered).}
  \item{reptech}{A vector of factors indicating the technical replicate to which each column belongs. Default is NULL (no experimental design is considered).}
  \item{nb.iter}{The number of iterations used for the multiple imputation method (see \code{\link{mi.mix}}).}
  \item{methodMCAR}{The method used for imputing MCAR data. If \code{methodMCAR="mle"} (default), then the \code{\link{impute.mle}} function is used (imputation using an EM algorithm). If \code{methodMCAR="pca"}, then the \code{\link{impute.PCA}} function is used (imputation using Principal Component Analysis). If \code{methodMCAR="rf"}, then the \code{\link{impute.RF}} function is used (imputation using Random Forest). Else, the \code{\link{impute.slsa}} function is used (imputation using Least Squares on nearest neighbours).}
  \item{methodMNAR}{The method used for imputing MNAR data. If \code{methodMNAR="igcda"} (default), then the \code{\link{impute.igcda}} function is used. Else, the \code{\link{impute.pa}} function is used.}
  \item{nknn}{The number of nearest neighbours used in the SLSA algorithm (see \code{\link{impute.slsa}}).}
  \item{selec}{A parameter to select a part of the dataset to find nearest neighbours between rows. This can be useful for big data sets (see \code{\link{impute.slsa}}).}
    \item{siz}{A parameter to select a part of the dataset to perform imputations with the MCAR-devoted algorithm. This can be useful for big data sets (see \code{\link{mi.mix}}).}
  \item{weight}{The way of weighting in the algorithm (see \code{\link{impute.slsa}}).}
  \item{ind.comp}{If \code{ind.comp=1}, only nearest neighbours without missing values are selected to fit linear models (see \code{\link{impute.slsa}}). Else, they can contain missing values.}
  \item{progress.bar}{If \code{TRUE}, a progress bar is displayed.}
  \item{x.step.mod}{The number of points in the intervals used for estimating the cumulative distribution functions of the mixing model in each column (see \code{\link{estim.mix}}).}
  \item{x.step.pi}{The number of points in the intervals used for estimating the proportion of MCAR values in each column (see \code{\link{estim.mix}}).}
  \item{nb.rei}{The number of initializations of the minimization algorithm used to estimate the proportion of MCAR values (see Details) (see \code{\link{estim.mix}}).}
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
First, a mixture model of MCAR and MNAR values is estimated in each column of \code{tab}. This model is used to estimate probabilities that each missing value is MCAR. Then, these probabilities are used to perform a multiple imputation strategy (see \code{\link{mi.mix}}). Rows with no value in a condition are imputed using the \code{\link{impute.pa}} function. More details and explanations can be bound in Giai Gianetto (2020).
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
res.sim=sim.data(nb.pept=2000,nb.miss=600,nb.cond=2);

#Imputation of the dataset noting the conditions to which the samples belong.
result=impute.mi(tab=res.sim$dat.obs, conditions=res.sim$conditions);

#Imputation of the dataset noting the conditions to which the samples belong
#and also their biological replicate, and using the SLSA method for the MCAR values
result=impute.mi(tab=res.sim$dat.obs, conditions=res.sim$conditions,
repbio=res.sim$repbio, methodMCAR = "slsa");

#For large data sets, the SLSA imputation can be accelerated thanks to the selec parameter
#and the siz parameter (see impute.slsa and mi.mix)
#but it may result in a less accurate data imputation. Note that selec has to be greater than siz.
#Here, nb.iter is fixed to 3
result1=impute.mi(tab=res.sim$dat.obs, conditions=res.sim$conditions, progress.bar=TRUE,
selec=400, siz=300, nb.iter=3, methodMCAR = "slsa", methodMNAR = "igcda");
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Imputation methods}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
