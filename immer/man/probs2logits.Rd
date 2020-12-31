%% File Name: probs2logits.Rd
%% File Version: 0.10

\name{probs2logits}
\alias{probs2logits}
\alias{logits2probs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Conversion of Probabilities into Logits
}
\description{
Converts probabilities into logits
}
\usage{
probs2logits(probs)

logits2probs(y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{probs}{
Vector containing probabilities
}
\item{y}{
Vector containing logits
}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{
A vector with logits or probabilities
}

%\references{
%% ~put references to the literature/web site here ~
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
# EXAMPLE 1: Probability-logit-conversions: a toy example
#############################################################################

# define vector of probabilities
probs <- c( .3, .25, .25, .2)
sum(probs)

# convert probabilities into logits
y <- immer::probs2logits( probs )
# retransform logits into probabilities
immer::logits2probs(y)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
