\name{predict.DStree}
\alias{predict.DStree}
\title{Predictions from a fitted DStree Object}
\usage{
\method{predict}{DStree}(object, data, ...)
}
\arguments{
  \item{object}{a fitted object of class \code{"DStree"}}

  \item{data}{data frame containing the values at which
  predictions are required. The predictors referred to in
  the right side of \code{formula(object)} as well as the
  'status' variable must be present by name in \code{data}. If
  missing, the fitted values are returned.}

  \item{...}{further arguments passed to or from other
  methods.}
}
\value{
A named list with the following elements:
The vector \code{MedSurv} contains the predicted median survival
times of the observations in \code{data}. The matrices \code{Surv} and
\code{Haz} contain the predicted probabilities and hazard rates, respectively,
where each row refers to an individual and each column refers to an observed time point (1,2,...).
}
\description{
Returns the predicted median survival time, hazard and
survival probabilities from a fitted \code{DStree} object.
}
\examples{
data(cost)
# Convert observed days to years
d.cost<-dis.cost(cost)

# Train Data
Train<-d.cost[1:300,]
# Test Data
Test<-d.cost[301:518,]

fit <- DStree(time~prevStroke+age+sex+alcohol+smoke,status="status",data=d.cost)

#Predictions from Test Data
predict(fit,Test)
}
\references{
Bou-Hamad I., Larocque D.,Ben-Ameur H., Masse L. C., Vitaro
F. and Tremblay R. E. (2009), Discrete-Time Survival Trees.
\emph{Canadian Journal of Statistics 37 (1), 17-32}.
}