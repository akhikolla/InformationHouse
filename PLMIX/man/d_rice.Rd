\name{d_rice}
\alias{d_rice}
\docType{data}
\title{Rice Voting Data (partial orderings)}
\description{The Rice Voting dataset (\code{d_rice}) collects the results of the 1992 election of a faculty member to serve on the Presidential Search Committee in the Rice University.
A total of \eqn{N=300} people casted their vote in the ballots by ranking the \eqn{K=5} candidates in the short list in a preferential manner. The dataset is composed of partial top orderings of varying lengths. Missing positions are denoted with zero entries.}
\usage{data(d_rice)}
\format{
Object of S3 class \code{c("top_ordering","matrix")} gathering a matrix of partial orderings with \eqn{N=300} rows and \eqn{K=5} columns. Each row lists the faculty members from the most-liked (\code{Rank_1}) to the least-liked (\code{Rank_5}) in a given voting ballot. 
}
\references{
Marcus, P., Heiser, W. J. and  D'Ambrosio, A. (2013). Comparison of heterogeneous probability models for ranking data, Master Thesis, Leiden University.

Baggerly, K. A. (1995). Visual estimation of structure in ranked data, PhD thesis, Rice University.

}
\examples{
data(d_rice)
head(d_rice)

## Subset of complete sequences
d_rice_compl=d_rice[rowSums(d_rice!=0)>=(ncol(d_rice)-1),]
head(d_rice_compl)

}
\keyword{datasets}

