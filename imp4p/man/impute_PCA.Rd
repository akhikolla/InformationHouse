\name{impute.PCA}
\alias{impute.PCA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Imputing missing values using Principal Components Analysis.
}
\description{
Imputing missing values using the algorithm proposed by Josse and Husson (2013). The function is based on the imputePCA function of the R package missMDA.
}
\usage{
impute.PCA(tab, conditions, ncp.max=5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tab}{A data matrix containing numeric and missing values. Each column of this matrix is assumed to correspond to an experimental sample, and each row to an identified peptide.}
  \item{conditions}{A vector of factors indicating the biological condition to which each sample belongs.}
  \item{ncp.max}{integer corresponding to the maximum number of components to test (used in the \code{estim_ncpPCA} function of R package missMDA).}
}
\details{
See Josse and Husson (2013) for the theory. It is built from functions proposed in the R package missMDA.
}
\value{
The input matrix \code{tab} with imputed values instead of missing values.
}
\references{
Josse, J & Husson, F. (2013). Handling missing values in exploratory multivariate data analysis methods. Journal de la SFdS. 153 (2), pp. 79-99.
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

#Imputation of missing values with PCA
dat.pca=impute.PCA(tab=res.sim$dat.obs,conditions=res.sim$condition);
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Imputation methods}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
