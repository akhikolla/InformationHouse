\name{plot.DStree}
\alias{plot.DStree}
\title{Plot a DStree Object}
\usage{
\method{plot}{DStree}(x, prob = "surv", select = NULL, ...)
}
\arguments{
  \item{x}{a fitted object of class \code{"DStree"}}

  \item{prob}{a string that indicates which probability should be
  plotted. \code{"surv"} draws survival probabilities,
  \code{"haz"} draws hazard probabilities.}

  \item{select}{a vector of strings that indicate terminal leaves
  for which the survival/hazard probabilities should be
  drawn. If not specified probabilities for every terminal leaf
  are drawn.}

  \item{...}{extra \code{\link{prp}} arguments to modify the
  appearance of the tree (from R package
  \bold{rpart.plot})}
}
\description{
This function plots a DStree object on the current graphic
device.  It visualizes the fitted tree as well as the
estimated survival/hazard probabilites for every terminal
node.
}
\examples{
data(cost)
# Discretize observed days to years
d.cost<-dis.cost(cost)

# Grow tree
fit <- DStree(time~prevStroke+age+sex+alcohol+smoke,status="status",data=d.cost)

#Plot tree and survival probabilities
plot(fit)
#survival probabilities are in the first plot

#Plot tree and hazard probabilities for terminal leaves 4 and 15
plot(fit,prob="haz",select=c("2","7"))
}

