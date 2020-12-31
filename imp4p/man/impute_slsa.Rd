\name{impute.slsa}
\alias{impute.slsa}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Imputing missing values using an adaptation of the LSimpute algorithm (Bo et al. (2004)) to experimental designs. This algorithm is named "Structured Least Squares Algorithm" (SLSA).
}
\description{
This function is an adaptation of the LSimpute algorithm (Bo et al. (2004)) to experimental designs usually met in MS-based quantitative proteomics.
}
\usage{
impute.slsa(tab, conditions, repbio=NULL, reptech=NULL, nknn=30, selec="all", weight="o",
ind.comp=1, progress.bar=TRUE)
  }
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tab}{A data matrix containing numeric and missing values. Each column of this matrix is assumed to correspond to an experimental sample, and each row to an identified peptide.}
  \item{conditions}{A vector of factors indicating the biological condition to which each sample belongs.}
  \item{repbio}{A vector of factors indicating the biological replicate to which each sample belongs. Default is NULL (no experimental design is considered).}
  \item{reptech}{A vector of factors indicating the technical replicate to which each sample belongs. Default is NULL (no experimental design is considered).}
  \item{nknn}{The number of nearest neighbours used in the algorithm (see Details).}
  \item{selec}{A parameter to select a part of the dataset to find nearest neighbours between rows. This can be useful for big data sets (see Details).}
  \item{weight}{The way of weighting in the algorithm (see Details).}
  \item{ind.comp}{If \code{ind.comp=1}, only nearest neighbours without missing values are selected to fit linear models (see Details). Else, they can contain missing values.}
  \item{progress.bar}{If \code{TRUE}, a progress bar is displayed.}
%%     ~~Describe \code{x} here~~
}
\details{
This function imputes the missing values condition by condition. The rows of the input matrix are imputed when they have at least one observed value in the considered condition. For the rows having only missing values in a condition, you can use the \code{\link{impute.pa}} function.

For each row, a similarity measure between the observed values of this row and the ones of the other rows is computed. The similarity measure which is used is the absolute pairwise correlation coefficient if at least three side-by-side values are observed, and the inverse of the euclidean distance between side-by-side observed values in the other cases.

For big data sets, this step can be time consuming and that is why the input parameter \code{selec} allows to select random rows in the data set. If \code{selec="all"}, then all the rows of the data set are considered; while if \code{selec} is a numeric value, for instance \code{selec=100}, then only 100 random rows are selected in the data set for computing similarity measures with each row containing missing values.

Once similarity measures are computed for a specific row, then the \code{nknn} rows with the highest similarity measures are considered to fit linear models and to predict several estimates for each missing value (see Bo et al. (2004)). If \code{ind.comp=1}, then only nearest neighbours without missing values in the condition are considered. However, unlike the original algorithm, our algorithm allows to consider the design of experiments that are specified in input through the vectors \code{conditions}, \code{repbio} and \code{reptech}. Note that \code{conditions} has to get a lower number of levels than \code{repbio}; and  \code{repbio} has to get a lower number of levels than \code{reptech}.

In the original algorithm, several predictions of each missing value are done from the estimated linear models and, then, they are weighted in function of their similarity measure and summed (see Bo et al. (2004)). In our algorithm, one can use the original weighting function of Bo et al. (2004) if \code{weight="o"}, i.e. \code{(sim^2/(1-sim^2+1e-06))^2} where \code{sim} is the similarity measure; or the weighting function \code{sim^weight} if \code{weight} is a numeric value.
}
\value{
The input matrix \code{tab} with imputed values instead of missing values.
}
\references{
Bo, T. H., Dysvik, B., & Jonassen, I. (2004). LSimpute: accurate estimation of missing values in microarray data with least squares methods. Nucleic acids research, 32(3), e34.
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

#Imputation of missing values with the slsa algorithm
dat.slsa=impute.slsa(tab=res.sim$dat.obs,conditions=res.sim$condition,repbio=res.sim$repbio);
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Imputation methods}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
