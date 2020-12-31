\name{impute.rand}
\alias{impute.rand}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Imputation of peptides with a random value.
}
\description{
For each row (peptide), this function imputes missing values by random values following a Gaussian distribution.
}
\usage{
impute.rand(tab, conditions)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tab}{A data matrix containing numeric and missing values. Each column of this matrix is assumed to correspond to an experimental sample, and each row to an identified peptide.}
  \item{conditions}{A vector of factors indicating the biological condition to which each column (experimental sample) belongs.}
  }
\details{
For each row (peptide), this function imputes missing values by random values following a Gaussian distribution centered on the mean of the observed values in the condition for the specific peptide and with a standard deviation equal to the first quartile of the distribution of the standard deviation the values observed for all the peptides. Rows with only missing values in a condition are not imputed (the \code{\link{impute.pa}} function can be used for this purpose).
}
\value{
The input matrix \code{tab} with imputed values instead of missing values.
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

#Imputation of the simulated data set with random values
data.rand=impute.rand(res.sim$dat.obs,res.sim$conditions);
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Imputation methods}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
