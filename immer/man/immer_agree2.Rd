%% File Name: immer_agree2.Rd
%% File Version: 0.18

\name{immer_agree2}
\alias{immer_agree2}
\alias{summary.immer_agree2}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Agreement Statistics for 2 Raters
}

\description{
Some agreement statistics for two raters, including raw
agreement, Scott's Pi, Cohen's Kappa, Gwets AC1 and
Aickens Alpha (see Gwet, 2010).
}
\usage{
immer_agree2(y, w=rep(1, nrow(y)), symmetrize=FALSE, tol=c(0, 1))

\method{summary}{immer_agree2}(object, digits=3,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
Data frame with responses for two raters
}
  \item{w}{
Optional vector of frequency weights
}
  \item{symmetrize}{
Logical indicating whether contingency table should
be symmetrized
}
  \item{tol}{
Vector of integers indicating tolerance for raw agreement
}
\item{object}{Object of class \code{immer_agree2}}
\item{digits}{Number of digits after decimal for rounding}
\item{\dots}{Further arguments to be passed}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
List with entries
\item{agree_raw}{Raw agreement}
\item{agree_stats}{Agreement statistics}
\item{agree_table}{Contingency table}
\item{marg}{Marginal frequencies}
\item{Pe}{Expected chance agreement probabilities}
\item{PH}{Probabilities for hard-to-classify
subjects according to Aicken}
\item{nobs}{Number of observations}
}

\references{
Gwet, K. L. (2010). \emph{Handbook of inter-rater reliability}.
Gaithersburg: Advanced Analytics.
}


%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
For more inter-rater agreement statistics see the \R packages
\pkg{agRee},
\pkg{Agreement},
\pkg{agrmt},
\pkg{irr},
\pkg{obs.agree},
\pkg{rel}.
}


\examples{
#############################################################################
# EXAMPLE 1: Dataset in Schuster & Smith (2006)
#############################################################################

data(data.immer08)
dat <- data.immer08

y <- dat[,1:2]
w <- dat[,3]
# agreement statistics
res <- immer::immer_agree2( y=y, w=w )
summary(res)
# extract some output values
res$agree_stats
res$agree_raw
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
