\name{impute.RF}
\alias{impute.RF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Imputing missing values using Random Forest.
}
\description{
Imputing missing values using the algorithm proposed by Stekhoven and Buehlmann (2012). The function is based on the missForest function of the R package missForest.
}
\usage{
impute.RF(tab, conditions,
                    maxiter = 10, ntree = 100, variablewise = FALSE,
                    decreasing = FALSE, verbose = FALSE,
                    mtry = floor(sqrt(ncol(tab))), replace = TRUE,
                    classwt = NULL, cutoff = NULL, strata = NULL,
                    sampsize = NULL, nodesize = NULL, maxnodes = NULL,
                    xtrue = NA, parallelize = c('no', 'variables', 'forests'))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tab}{A data matrix containing numeric and missing values. Each column of this matrix is assumed to correspond to an experimental sample, and each row to an identified peptide.}
  \item{conditions}{A vector of factors indicating the biological condition to which each sample belongs.}
  \item{maxiter}{parameter of the \code{missForest} function (missForest R package).}
  \item{ntree}{parameter of the \code{missForest} function (missForest R package).}
  \item{variablewise}{parameter of the \code{missForest} function (missForest R package).}
  \item{decreasing}{parameter of the \code{missForest} function (missForest R package).}
  \item{verbose}{parameter of the \code{missForest} function (missForest R package).}
  \item{mtry}{parameter of the \code{missForest} function (missForest R package).}
  \item{replace}{parameter of the \code{missForest} function (missForest R package).}
  \item{classwt}{parameter of the \code{missForest} function (missForest R package).}
  \item{cutoff}{parameter of the \code{missForest} function (missForest R package).}
  \item{strata}{parameter of the \code{missForest} function (missForest R package).}
  \item{sampsize}{parameter of the \code{missForest} function (missForest R package).}
  \item{nodesize}{parameter of the \code{missForest} function (missForest R package).}
  \item{maxnodes}{parameter of the \code{missForest} function (missForest R package).}
  \item{xtrue}{parameter of the \code{missForest} function (missForest R package).}
  \item{parallelize}{parameter of the \code{missForest} function (missForest R package).}
}
\details{
See Stekhoven and Buehlmann (2012) for the theory. It is built from functions proposed in the R package missForest.
}
\value{
The input matrix \code{tab} with imputed values instead of missing values.
}
\references{
Stekhoven, D.J. and Buehlmann, P. (2012), 'MissForest - nonparametric missing value imputation for mixed-type data', Bioinformatics, 28(1) 2012, 112-118, doi: 10.1093/bioinformatics/btr597
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

#Imputation of missing values with Random Forest
dat.rf=impute.RF(tab=res.sim$dat.obs,conditions=res.sim$condition);
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Imputation methods}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
