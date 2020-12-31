\name{fast_apply_nb_na}
\alias{fast_apply_nb_na}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function similar to the function \code{apply(X,dim,function(x)sum(is.na(x)))}.
}
\description{
This function is similar to the function \code{apply(X,dim,function(x)sum(is.na(x)))} but written thanks to the Rcpp package and therefore faster than \code{apply(X,dim,function(x)sum(is.na(x)))}.
}
\usage{
fast_apply_nb_na(X, dim)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
A data matrix containing numeric and missing values.
}
  \item{dim}{
A numeric value: 1 if the number of missing values has to be computed for each row of \code{X}, or 2 if it has to be computed for each column of \code{X}.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A numeric vector containing the number of missing values in either each row or each column of \code{X}.
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
##    .Call("imp4p_fast_apply_nb_na", PACKAGE = "imp4p", X, dim)
##  }
##
## You can compare the execution time with a traditional apply function by
## library(rbenchmark)
## res.sim=sim.data(nb.pept=2000,nb.miss=600);
## benchmark(fast_apply_nb_na(res.sim$dat.obs, 1),
##           apply(res.sim$dat.obs,1,function(x)sum(is.na(x))))
##
##
##

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Cpp wrapper function}% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
