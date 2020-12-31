\name{d_carconf}
\alias{d_carconf}
\docType{data}
\title{Car Configurator Data (partial orderings)}
\description{The Car Configurator dataset (\code{d_carconf}) came up from a marketing study aimed at investigating customer preferences toward different car features.
A sample of \eqn{N=435} customers were asked to construct their car by using an online configurator system and choose among \eqn{K=6} car modules in order of preference. The car features are labeled as:
 1 = price, 2 = exterior design, 3 = brand, 4 = technical equipment, 5 = producing country and 6 = interior design.
The survey did not require a complete ranking elicitation, therefore the dataset is composed of partial top orderings of varying lengths. Missing positions are denoted with zero entries.}
\usage{data(d_carconf)}
\format{
Object of S3 class \code{c("top_ordering","matrix")} gathering a matrix of partial orderings with \eqn{N=435} rows and \eqn{K=6} columns. Each row lists the car features from the most important (\code{Rank_1}) to the least important (\code{Rank_6}) for a given customer.
}
\references{
Mollica, C. and Tardella, L. (2017). Bayesian Plackett-Luce mixture models for partially ranked data. \emph{Psychometrika}, \bold{82}(2), pages 442--458, ISSN: 0033-3123, DOI: 10.1007/s11336-016-9530-0.

Hatzinger, R. and Dittrich, R. (2012). Prefmod: An R package for modeling preferences based on paired comparisons, rankings, or ratings. \emph{Journal of Statistical Software}, \bold{48}(10), pages 1--31.

Dabic, M. and Hatzinger, R. (2009). Zielgruppenadaequate Ablaeufe in Konfigurationssystemen - eine empirische Studie im Automobilmarkt - Partial Rankings. In Hatzinger, R., Dittrich, R. and  Salzberger, T. (eds), \emph{Praeferenzanalyse mit R: Anwendungen aus Marketing, Behavioural
Finance und Human Resource Management}. Wien: Facultas.

}
\examples{
data(d_carconf)
head(d_carconf)

## Subset of complete sequences
d_carconf_compl=d_carconf[rowSums(d_carconf!=0)>=(ncol(d_carconf)-1),]
head(d_carconf_compl)
}
\keyword{datasets}
