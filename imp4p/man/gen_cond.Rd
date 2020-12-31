\name{gen.cond}
\alias{gen.cond}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function allowing to create a vector indicating the membership of each sample to a condition.
}
\description{
This function creates a vector of factors where each element refers to a condition to which a sample belongs.}
\usage{
gen.cond(nb_cond=2,nb_sample=5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nb_cond}{Number of biological conditions.}
  \item{nb_sample}{Number of samples in each condition.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A vector of factors of length nb_cond*nb_sample.
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

cond=gen.cond(nb_cond=2,nb_sample=6)
#[1] 1 1 1 1 1 1 2 2 2 2 2 2
#Levels: 1 2
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Simulated data}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
