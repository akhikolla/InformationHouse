\name{d_gaming}
\alias{d_gaming}
\docType{data}
\title{Gaming Platforms Data (complete orderings)}
\description{The Gaming Platforms dataset (\code{d_gaming}) collects the results of a survey conducted on a sample of \eqn{N=91} Dutch students, who were asked to rank \eqn{K=6} gaming platforms in order of preference, namely: 1 = X-Box, 2 = PlayStation, 3 = PSPortable, 4 = GameCube, 5 = GameBoy and 6 = Personal Computer. The dataset is composed of complete orderings. 
}
\usage{data(d_gaming)}
\format{
Object of S3 class \code{c("top_ordering","matrix")} gathering a matrix of complete orderings with \eqn{N=91} rows and \eqn{K=6} columns. Each row lists the gaming platforms from the most-liked (\code{Rank_1}) to the least-liked (\code{Rank_6}) for a given student.
}
\references{
Fok, D., Paap, R. and Van Dijk, B. (2012). A Rank-Ordered Logit Model With Unobserved Heterogeneity In Ranking Capatibilities. \emph{Journal of Applied Econometrics}, \bold{27}(5), pages 831--846.
}
\source{
The Gaming Platforms dataset in .csv format can be downloaded from \url{http://qed.econ.queensu.ca/jae/2012-v27.5/fok-paap-van_dijk/}. The .csv files contains the preference data in ranking format and some covariates collected for each student.}
\examples{
data(d_gaming)
head(d_gaming)
}
\keyword{datasets}

