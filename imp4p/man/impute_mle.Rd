\name{impute.mle}
\alias{impute.mle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Imputing missing values using a maximum likelihood estimation (MLE).
}
\description{
Imputing missing values using the EM algorithm proposed in section 5.4.1 of Schafer (1997). The function is based on the imp.norm function of the R package norm.
}
\usage{
impute.mle(tab, conditions)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tab}{A data matrix containing numeric and missing values. Each column of this matrix is assumed to correspond to an experimental sample, and each row to an identified peptide.}
  \item{conditions}{A vector of factors indicating the biological condition to which each sample belongs.}
}
\details{
See section 5.4.1 of Schafer (1997) for the theory. It is built from functions proposed in the R package norm.
}
\value{
The input matrix \code{tab} with imputed values instead of missing values.
}
\references{
Schafer, J. L. (1997). Analysis of incomplete multivariate data. Chapman and Hall/CRC.
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

#Imputation of missing values with the mle algorithm
dat.mle=impute.mle(tab=res.sim$dat.obs,conditions=res.sim$condition);
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Imputation methods}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
