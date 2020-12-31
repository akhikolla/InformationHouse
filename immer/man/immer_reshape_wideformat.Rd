%% File Name: immer_reshape_wideformat.Rd
%% File Version: 0.18

\name{immer_reshape_wideformat}
\alias{immer_reshape_wideformat}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Creating a Rating Dataset in Wide Format
}
\description{
Converts a rating dataset from a long format into a wide format.
}
\usage{
immer_reshape_wideformat(y, pid, rater, Nmin_ratings=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
Vector or a data frame containing ratings
}
  \item{pid}{
Person identifier
}
  \item{rater}{
Rater identifier
}
  \item{Nmin_ratings}{
Minimum number of ratings used for selection
}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
Data frame with ratings. Each row corresponds to a person, and each
of the columns (except the first one containing the person
identifier) to one rater.
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
# EXAMPLE 1: Reshaping ratings of one variable into wide format
#############################################################################

data(data.immer03)
dat <- data.immer03

# select variable "b" and persons which have at least 10 ratings
dfr <- immer::immer_reshape_wideformat( y=dat$b2, pid=dat$idstud, rater=dat$rater,
                    Nmin_ratings=10 )
head(dfr)

#############################################################################
# EXAMPLE 2: Reshaping ratings of a data frame into wide format
#############################################################################

data(data.immer07)
dat <- data.immer07

#*** Dataset 1: Wide format for item I1
dfr1 <- immer::immer_reshape_wideformat( dat$I1, rater=dat$rater, pid=dat$pid)

#*** Dataset 2: Wide format for four items I1, I2, I3 and I4
dfr2 <- immer::immer_reshape_wideformat( dat[, paste0("I",1:4) ],
              rater=dat$rater, pid=dat$pid )
str(dfr2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Utility function}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
