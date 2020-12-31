\name{simMLgvar}
\alias{simMLgvar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate graphical VAR data of multiple subjects
}
\description{
See arxiv.org/abs/1609.04156 for details.
}
\usage{
simMLgvar(nTime, nVar, nPerson, propPositive = 0.5, kappaRange = c(0.25, 0.5), 
          betaRange = c(0.25, 0.5), betweenRange = c(0.25, 0.5), 
          rewireWithin = 0, betweenVar = 1, withinVar = 0.25, 
          temporalOffset = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nTime}{
Number of time points per subject
}
  \item{nVar}{
Number of variables
}
  \item{nPerson}{
Number of subjects
}
  \item{propPositive}{
Proportion of positive edges
}
  \item{kappaRange}{
Range of partial contemporaneous correlation coefficients
}
  \item{betaRange}{
Range of temporal coefficients
}
  \item{betweenRange}{
Range of partial between-subjects coefficients
}
  \item{rewireWithin}{
Rewiring probability of contemporaneous networks
}
  \item{betweenVar}{
Between-subjects variabce
}
  \item{withinVar}{
Contemporaneous variance
}
  \item{temporalOffset}{
Specifies the temporal network. Setting this to 2 connects X_i to X_(i+2)
}
}
\value{
A \code{"simMLgvar"} object with the following elements:
\item{data }{Generated dataset}
\item{fixedKappa }{Fixed inverse contemporaneous covariance matrix}
\item{fixedPCC }{Fixed contemporaneous partial correlation network}
\item{fixedBeta }{Fixed temporal network}
\item{fixedPDC }{Fixed standardized temporal network}
\item{between }{Fixed between-subjects network}
\item{means }{True means}
\item{personData }{Dataset split per person}
\item{idvar }{String indicating the id variable}
\item{vars }{Vector of strings indicating the variables}
}
\author{
Sacha Epskamp <mail@sachaepskamp.com>
}