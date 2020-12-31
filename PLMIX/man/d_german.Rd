\name{d_german}
\alias{d_german}
\docType{data}
\title{German Sample Data (complete orderings)}
\description{The German Sample dataset (\code{d_german}) is part of a comparative cross-sectional study on political actions and mass participation involving five Western countries. The dataset regards a sample of \eqn{N=2262} German respondents who were asked to rank \eqn{K=4} political goals in order of desirability, namely: 1 = maintaining order in the nation, 2 = giving people more say in the decisions of government, 3 = fighting rising prices and 4 = protecting freedom of speech. The dataset is composed of complete orderings. 
}
\usage{data(d_german)}
\format{
Object of S3 class \code{c("top_ordering","matrix")} gathering a matrix of complete orderings with \eqn{N=2262} rows and \eqn{K=4} columns. Each row lists the political goals from the most desiderable (\code{Rank_1}) to the least desiderable (\code{Rank_4}) for a given respondent.
}
\references{
Croon, M. A. (1989). Latent class models for the analysis of rankings. In De Soete, G., Feger, H. and Klauer, K. C. (eds), \emph{New Developments in Psychological Choice Modeling}, pages 99--121. North-Holland: Amsterdam.

Barnes, S. H. et al. (1979). Political action. Mass participation in five Western democracies. London: Sage.
}
\examples{
data(d_german)
head(d_german)
}
\keyword{datasets}
