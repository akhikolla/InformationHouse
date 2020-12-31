\name{predict.robmixglm}
\alias{predict.robmixglm}
\title{Predict Method for robmixglm}
\usage{
\method{predict}{robmixglm}(object,  newdata  =  NULL, 
            type  =  c("link",  "response"),  \ldots)
}
\arguments{
  \item{object}{a fitted object of class inheriting from \code{robmixglm}.}
  \item{newdata}{optionally,  a data frame in which to look for variables with
    which to predict.  If omitted,  the fitted linear predictors are used.}
  \item{type}{the type of prediction required.  The default \code{link} is on the
    scale of the linear predictors,  while the alternative \code{response}
    is on the scale of the response variable.}  
  \item{\ldots}{Other parameters. (not used)}
}
\description{
  Obtains predictions from a fitted robust mixture generalized linear model object.
}
\details{
  If \code{newdata} is omitted the predictions are based on the data
  used for the fit.  In that case how cases with missing values in the
  original fit is determined by the \code{na.action} argument of that
  fit.  If \code{na.action  =  na.omit} omitted cases will not appear in
  the residuals,  whereas if \code{na.action  =  na.exclude} they will
  appear (in predictions and standard errors),  with residual value
  \code{NA}.  See also \code{\link{napredict}}.
}
\value{
  A vector predicted linear predictors or response. For \code{binomial} the
  resonse is the predicted proportion.
  
}

\examples{
\donttest{
library(MASS)
data(forbes)
forbes.robustmix <- robmixglm(100*log10(pres)~bp, data = forbes, cores = 1)
plot(forbes$bp, forbes$pres)
preddata <- data.frame(bp = seq(from  =  min(forbes$bp), to  =  max(forbes$bp),  by  =  0.01))
# convert to original scale
preddata$predpres <-10^(predict(forbes.robustmix, newdata = preddata)/100)
lines(preddata$bp, preddata$predpres, col = "red")
}}
\keyword{models}
