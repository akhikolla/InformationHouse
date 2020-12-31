\name{mlGraphicalVAR}
\alias{mlGraphicalVAR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Pooled and individual graphical VAR estimation
}
\description{
This function fits fixed effect temporal and contemporaneous networks over multiple subjects and runs separate graphical VAR models per subject. The algorithm does: (1) pool all data, within-subject center variables and run \code{\link{graphicalVAR}} to obtain fixed effects, (2) run \code{\link[qgraph]{EBICglasso}} on subject means to obtain a between-subjects network, (3) run \code{\link{graphicalVAR}} on data of every subject to obtain individual networks.  See arxiv.org/abs/1609.04156 for more details.
}
\usage{
mlGraphicalVAR(data, vars, beepvar, dayvar, idvar, scale = TRUE, 
              centerWithin = TRUE, gamma = 0.5, verbose = TRUE, 
              subjectNetworks = TRUE, lambda_min_kappa_fixed = 0.001, 
              lambda_min_beta_fixed = 0.001, lambda_min_kappa = 0.05, 
              lambda_min_beta = lambda_min_kappa, lambda_min_glasso = 0.01, 
              ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Data frame
}
  \item{vars}{
Vectors of variables to include in the analysis
}
  \item{beepvar}{
String indicating assessment beep per day (if missing, is added). Adding this argument will cause non-consecutive beeps to be treated as missing!
}
  \item{dayvar}{
String indicating assessment day. Adding this argument makes sure that the first measurement of a day is not regressed on the last measurement of the previous day. IMPORTANT: only add this if the data has multiple observations per day.
}
  \item{idvar}{
String indicating the subject ID
}
\item{scale}{Logical, should variables be standardized before estimation?}
  \item{centerWithin}{
Logical, should subject data be within-person centered before estimating fixed effects?
}
  \item{gamma}{
EBIC tuning parameter.
}
  \item{verbose}{
Logical indicating if console messages and the progress bar should be shown.
}
  \item{subjectNetworks}{
\code{TRUE} to estimate all subject numbers, or a vector with IDs of which subject numbers should be estimated.
}
  \item{lambda_min_kappa_fixed}{
Multiplier of maximal tuning parameter
}
  \item{lambda_min_beta_fixed}{
Multiplier of maximal tuning parameter
}
  \item{lambda_min_kappa}{
Multiplier of maximal tuning parameter
}
  \item{lambda_min_beta}{
Multiplier of maximal tuning parameter
}
  \item{lambda_min_glasso}{
Multiplier of maximal tuning parameter
}
  \item{\dots}{
Arguments sent to \code{\link{graphicalVAR}}
}
}
\value{
A \code{"mlGraphicalVAR"} object with the following elements:
\item{fixedPCC }{Estimated fixed effects (partial contemporaneous correlations) of contemporaneous effects}
\item{fixedPDC }{Estimated fixed effects (partial directed correlations) of temporal effects}
\item{fixedResults }{Full object of pooled data estimation (fixed effects)}
\item{betweenNet }{Estimated between-subjects network (partial correlations)}
\item{ids }{Vector of subject IDs}
\item{subjectPCC }{List of estimated individual contemporaneous networks}
\item{subjectPDC }{List of estimated individual directed networks}
\item{subjecResults }{List of full results of individual estimations}
}
\references{
Epskamp, S., Waldorp, L. J., M\~{o}ttus, R., & Borsboom, D. Discovering Psychological Dynamics: The Gaussian Graphical Model in Cross-sectional and Time-series Data.
}
\author{
Sacha Epskamp <mail@sachaepskamp.com>
}
\seealso{
\code{\link{graphicalVAR}}
}
\examples{
\dontrun{
# Simulate data:
Sim <- simMLgvar(nTime = 50, nPerson = 20, nVar = 3)

# Estimate model:
Res <- mlGraphicalVAR(Sim$data, vars = Sim$vars, idvar = Sim$idvar)

layout(t(1:2))
library("qgraph")

# Temporal fixed effects
qgraph(Res$fixedPDC, title = "Estimated fixed PDC", layout = "circle")
qgraph(Sim$fixedPDC, title = "Simulated fixed PDC", layout = "circle")

# Contemporaneous fixed effects
qgraph(Res$fixedPCC, title = "Estimated fixed PCC", layout = "circle")
qgraph(Sim$fixedPCC, title = "Simulated fixed PCC", layout = "circle")
}
}
