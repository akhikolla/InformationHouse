\name{impute.pa}
\alias{impute.pa}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Imputation of peptides having no value in a biological condition (present in a condition / absent in another).
}
\description{
This function imputes missing values by small values.
}
\usage{
impute.pa(tab, conditions, q.min = 0.025, q.norm = 3, eps = 0,
distribution = "unif", param1 = 3, param2 = 1, R.q.min=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tab}{A data matrix containing numeric and missing values. Each column of this matrix is assumed to correspond to an experimental sample, and each row to an identified peptide.}
  \item{conditions}{A vector of factors indicating the biological condition to which each column (experimental sample) belongs.}
  \item{q.min}{A quantile value of the observed values allowing defining the maximal value which can be generated. This maximal value is defined by the quantile \code{q.min} of the observed values distribution minus \code{eps}. Default is 0.025 (the maximal value is the 2.5 percentile of observed values minus \code{eps}).}
  \item{q.norm}{A quantile value of a normal distribution allowing defining the minimal value which can be generated. Default is 3 (the minimal value is the maximal value minus qn*median(sd(observed values)) where sd is the standard deviation of a row in a condition).}
  \item{eps}{A value allowing defining the maximal value which can be generated. This maximal value is defined by the quantile \code{q.min} of the observed values distribution minus \code{eps}. Default is 0.}
  \item{distribution}{Distribution used to generated missing values. You have the choice between "unif" for the uniform distribution, "beta" for the Beta distribution or "dirac" for the Dirac distribution. Default is "unif".}
  \item{param1}{Parameter \code{shape1} of the Beta distribution.}
  \item{param2}{Parameter \code{shape2} of the Beta distribution.}
  \item{R.q.min}{Parameter used for the Dirac distribution. In this case, all the missing values are imputed by a single value which is equal to \code{R.q.min*quantile(tab[,j], probs=q.min, na.rm=T)}. Default is 1 : the imputed value is the \code{qmin} quantile of observed values.}
}
\details{
This function replaces the missing values in a column by random draws from a specified distribution. The value of \code{eps} can be interpreted as a minimal fold-change value above which the present/absent peptides appear.
}
\value{
A list composed of :

- \code{tab.imp} : the input matrix \code{tab} with imputed values instead of missing values.

- \code{para} : the parameters of the distribution which has been used to impute.
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

#Imputation of the simulated data set with small values
data.small.val=impute.pa(res.sim$dat.obs,res.sim$conditions);
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Imputation methods}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
