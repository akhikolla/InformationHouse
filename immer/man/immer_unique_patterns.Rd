%% File Name: immer_unique_patterns.Rd
%% File Version: 0.09

\name{immer_unique_patterns}
\alias{immer_unique_patterns}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extracts Unique Item Response Patterns
}
\description{
Extracts unique item response patterns.
}
\usage{
immer_unique_patterns(dat, w=rep(1, nrow(dat)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
Data frame containing integer item responses
}
  \item{w}{
Optional vector of weights
}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{
A list with entries
  \item{y}{Data frame with unique item response patterns}
  \item{w}{Vector of frequency weights}
  \item{y_string}{Item response pattern coded as a string}
}
%\references{
%% ~put references to the literature/web site here ~
%}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See \code{\link[mirt:expand.table]{mirt::expand.table}} for back-converting
unique item response patterns into a data frame with item responses.
}

\examples{
#############################################################################
# EXAMPLE 1: Unique item response patterns data.read
#############################################################################

data( data.read, package="sirt")
dat <- data.read

# extract item response patterns
res <- immer::immer_unique_patterns( dat)

\dontrun{
# back-conversion with expand.table
dat2 <- mirt::expand.table( cbind( res$y,  res$w ) )
# check correctness
colMeans(dat)
colMeans(dat2)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
