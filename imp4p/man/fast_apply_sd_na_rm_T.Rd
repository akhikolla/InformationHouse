\name{fast_apply_sd_na_rm_T}
\alias{fast_apply_sd_na_rm_T}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function similar to the function \code{apply(X,dim,sd,na.rm=TRUE)}.
}
\description{
This function is similar to the function \code{apply(X,dim,sd,na.rm=TRUE)} but written thanks to the Rcpp package and therefore faster than \code{apply(X,dim,sd,na.rm=TRUE)}.
}
\usage{
fast_apply_sd_na_rm_T(X, dim)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
A data matrix containing numeric and missing values.
}
  \item{dim}{
A numeric value: 1 if the standard deviation has to be computed for each row of \code{X}, or 2 if the standard deviation has to be computed for each column of \code{X}.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A numeric vector containing the standard deviation of observed values in either each row or each column of \code{X}.
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
%\code{\link{impute.slsa}}
%}
\examples{
## The function is currently defined as
##function (X, dim)
##{
##    .Call("imp4p_fast_apply_sd_na_rm_T", PACKAGE = "imp4p", X,
##        dim)
##  }
##
## You can compare the execution time with a traditional apply function by
## library(rbenchmark)
## res.sim=sim.data(nb.pept=2000,nb.miss=600);
## benchmark(fast_apply_sd_na_rm_T(res.sim$dat.obs, 1),
##           apply(res.sim$dat.obs,1,sd,na.rm=TRUE))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Cpp wrapper function}% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
