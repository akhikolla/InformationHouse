\name{d_occup}
\alias{d_occup}
\docType{data}
\title{Occupation Data (complete orderings)}
\description{The Occupation dataset (\code{d_occup}) came up from a survey conducted on graduates from the Technion-Insrael Institute of Tecnology.
A sample of \eqn{N=143} graduates were asked to rank \eqn{K=10} professions according to the perceived prestige. The occupations are labeled as:
 1 = faculty member, 2 = owner of a business, 3 = applied scientist, 4 = operations researcher, 5 = industrial engineer, 6 = manager, 7 = mechanical engineer, 8 = supervisor, 9 = technician and 10 = foreman. The dataset is composed of complete orderings.}
\usage{data(d_occup)}
\format{
Object of S3 class \code{c("top_ordering","matrix")} gathering a matrix of complete orderings with \eqn{N=143} rows and \eqn{K=10} columns. Each row lists the professions from the most-liked (\code{Rank_1}) to the least-liked (\code{Rank_10}) for a given graduate.
}
\references{

Cohen, A. and Mallows, C. L. (1983). Assessing goodness of fit of ranking models to data. \emph{Journal of the Royal Statistical Society: Series D (The Statistician)}, \bold{32}(4), pages 361--374, ISSN: 0039-0526.

Cohen, A. (1982). Analysis of large sets of ranking data. \emph{Communications in Statistics -- Theory and Methods}, \bold{11}(3), pages 235--256.

Goldberg, A. I. (1976). The relevance of cosmopolitan/local orientations to professional values and behavior. \emph{Sociology of Work and Occupations}, \bold{3}(3), pages 331--356.

}
\examples{
data(d_occup)
head(d_occup)
}
\keyword{datasets}
