\name{dis.cost}
\alias{dis.cost}
\title{Convert observed time points in Copenhagen Stroke Study data set}
\description{
This function rounds up the observed number of days in the \code{cost} data set (contained in R package \bold{pec}) to full years.
}
\usage{
dis.cost(cost)
}
\arguments{
\item{cost}{data set}}
\value{
A data frame that contains the same variables as the \code{cost} data frame in R package \bold{pec}. The only difference between \code{cost} and the output of \code{dis.cost} are the units of the event times. In \code{cost}, event times are measured in days whereas \code{dis.cost} results in a data frame with event times measured in years.}
\references{
Joergensen HS, Nakayama H, Reith J, Raaschou HO, and Olsen TS.
Acute stroke with atrial fibrillation. The Copenhagen Stroke Study. Stroke,
27(10):1765-9, 1996.

Mogensen UB, Ishwaran H, and Gerds TA. Evaluating random forests for
survival analysis using prediction error curves. Journal of Statistical Software,
50(11), 2012.
}
\examples{
data(cost)
d.cost <- dis.cost(cost)}
