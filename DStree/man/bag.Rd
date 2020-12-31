\name{bag}
\alias{bag}
\title{Bagging Discrete-Time Survival Trees}
\usage{
bag(formula, data, status, nBoot = 10, minbucket = 40)
}
\arguments{
  \item{formula}{a \link{formula} with a response but no
  interaction terms.  The response variable represents the
  observed survival times.}

  \item{status}{integer/string representing the column
  number/name of the 'status' variable in the data frame.
  'status' must be a binary integer variable indicating whether the event
  occured (\code{=1}) or whether the event was censored/did not
  occur (\code{=0}).}

  \item{data}{data frame that contains all variables stated
  in the \code{formula} argument, as well as the
  \code{status} variable.}

  \item{nBoot}{an integer referring to the number of bootstrap
  replications.}

  \item{minbucket}{the minimum number of observations in
  any terminal node.}
}
\value{
The \code{bag} function returns an object of class 'DStreebag' which is a list
containing the following objects:

The vector MedSurv contains the predicted, averaged median survival
times. The matrices Surv and Haz contain the
predicted, averaged probabilities and hazard rates, respectively, where each row refers to an individual and each column refers to an observed time point (1,2,...). The list element 'trees' contains the \code{nBoot} unpruned trees of class 'DStree'. The list elements 'minbucket' and 'nboot' have
the same meaning as above.
}
\description{
Function for bagging discrete-time survival trees.
}
\details{
The \code{bag} function grows \code{nBoot} unpruned trees from
bootstrap samples. For each fitted tree the median survival
time, as well as the predicted survival probabilities and hazard rates
of the individuals in \code{data} are returned.
}
\examples{
data(cost)
## Discretize observed days to years
d.cost <- dis.cost(cost)# Bagging Tree
pred <- bag(time~prevStroke+age+sex+alcohol+smoke,status="status",data=d.cost,nBoot=50)

# Predicted, averaged probabilities and median survival times for each individual
pred$MedSurv
pred$Surv
pred$Haz
}
\references{
Bou-Hamad I., Larocque D., Ben-Ameur H., Masse L. C., Vitaro
F. and Tremblay R. E. (2009), Discrete-Time Survival Trees.
\emph{Canadian Journal of Statistics 37 (1), 17-32}.

Hothorn T., Lausen B., Benner A. and Radespiel-Troeger M.
(2004), Bagging Survival Trees. \emph{Statistics in
Medicine 23 (1), 77-91}.
}

