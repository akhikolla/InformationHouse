%% File Name: mnlfa_expand_to_list.Rd
%% File Version: 0.04

\name{mnlfa_expand_to_list}
\alias{mnlfa_expand_to_list}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Expands Input Into a List
}

\description{
Expands an input into a list.
}
\usage{
mnlfa_expand_to_list(x, names_list)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An R object
}
  \item{names_list}{
Names of the list
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A list
}

%\references{
%% ~put references to the literature/web site here ~
%}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
#############################################################################
# EXAMPLE 1: Test example
#############################################################################

formula_int <- ~1+female+age
items <- paste0("I",1:12)
formula_int <- mnlfa::mnlfa_expand_to_list(x=formula_int, names_list=items)
formula_int[[1]] <- ~0
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
