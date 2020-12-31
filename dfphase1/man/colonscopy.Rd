\name{colonscopy}
\alias{colonscopy}
\docType{data}
\title{Colonscopy Times}
\description{
  This data set contains the colonscopy times (minutes)
  for 30 subgroups of 5 patients given in Allison Jones-Farmer \emph{et al.} (2009).
}
\usage{data(colonscopy)}
\format{
  A 5x30 matrix.
}
\references{
    L. A. Jones-Farmer, V. Jordan, C. W. Champs (2009)
    \dQuote{Distribution-free Phase I control charts for subgroup location},
    \emph{Journal of Quality Technology}, \bold{41}, pp. 304--316.
}
\examples{
data(colonscopy)
phase1Plot(colonscopy)
}
\keyword{datasets}
