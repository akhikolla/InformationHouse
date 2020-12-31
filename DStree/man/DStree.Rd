\name{DStree}
\alias{DStree}
\title{Fit a discrete-time survival tree}
\usage{
DStree(formula, status, data, control = control, weights = NULL)
}
\arguments{
  \item{formula}{a \code{\link{formula}} with response but no
  interaction terms. The response variable represents the
  observed survival times.}

  \item{status}{integer/string representing the column
  number/name of the 'status' variable in the data frame.
  'status' must be a binary integer variable indicating whether the event
  occured (\code{=1}) or whether the event was censored/did not
  occur (\code{=0}).}

  \item{data}{data frame that contains all variables stated
  in the \code{formula} argument, as well as the 'status'
  variable.}

  \item{control}{a list of options that control
  the specification of the \code{DStree} algorithm. See \code{\link{DStree.control}}.}

  \item{weights}{a vector of optional case weights.}
  
}
\value{
An object of class \code{DStree}. For details see
\code{\link{DStree.object}}.
}
\description{
This function builds decision trees for discrete, right-censored
survival data. The fitted tree estimates the hazard and survival probabilities for every terminal node, as well as the median survival time.
}
\examples{
##Build tree
fit<- DStree(spell~ui+age+tenure+logwage,status="censor1",data=UnempDur,control=list(cp=0))
plot(fit)
}
\references{
Bou-Hamad I., Larocque D., Ben-Ameur H., Masse L. C., Vitaro
F. and Tremblay R. E. (2009), Discrete-Time Survival Trees.
\emph{Canadian Journal of Statistics 37 (1), 17-32}.
}

