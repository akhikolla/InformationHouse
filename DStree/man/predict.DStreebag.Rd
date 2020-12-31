\name{predict.DStreebag}
\alias{predict.DStreebag}
\title{Prediction from bagged Discrete-Time Survival Trees}
\usage{
\method{predict}{DStreebag}(object, data, ...)
}
\arguments{
  \item{object}{a fitted object of class \code{DStreebag}}

  \item{data}{data frame containing the values at which
  predictions are required. The predictors referred to in
  the right side of \code{formula(object)} as well as the 'status'
  variable must be present by name in \code{newdata}.}

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
Predict the median survival time, hazard and survival
probabilities from a \code{DStreebag} object.
}
\examples{
data(cost)
## Discretize observed days to years
d.cost<-dis.cost(cost)

# Bagging Tree
bag <- bag(time~prevStroke+age+sex+alcohol+smoke,status="status",data=d.cost[1:330,],nBoot=10)
pred <- predict(bag,newdata=d.cost[331:518,])
pred$MedSurv
}
\references{
Bou-Hamad I., Larocque D.,Ben-Ameur H., Masse L. C., Vitaro
F. and Tremblay R. E. (2009), Discrete-Time Survival Trees.
\emph{Canadian Journal of Statistics 37 (1), 17-32}.

Hothorn T., Lausen B., Benner A. and Radespiel-Troeger M.
(2004), Bagging Survival Trees. \emph{Statistics in
Medicine 23 (1), 77-91}.
}

