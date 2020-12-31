\name{predict.augSIMEX}
\alias{predict.augSIMEX}

\title{
Predict Method for the model fits by augSIMEX
}
\description{
This function returns the predictions and optionally estimates the standard errors of the predictions from a fitted augSIMEX object.
}
\usage{
## S3 method for class 'augSIMEX'

\method{predict}{augSIMEX}(object, newdata = NULL, type = c("link", "response","terms"), 
se.fit = FALSE, dispersion = NULL, terms = NULL, na.action = na.pass, ...) 
}

\arguments{
  \item{object}{
the ``augSIMEX" object gotten from \code{\link[augSIMEX]{augSIMEX}} function.
}
  \item{newdata}{
 An optional data frame in which to look for variables with which to predict. If not specified, the prediction will be conducted on the original main data.
}
  \item{type}{
 the type of prediction needed. The default is on the scale of the linear predictors; the alternative "\code{response}" is on the scale of the response variable. The setting is similar to the \code{\link[stats]{predict.glm}}. The "\code{terms}" option returns a matrix including the fitted values of each term in the model formula on the linear predictor scale.
}
  \item{se.fit}{
a logical variable indicating if standard errors are required.
}
  \item{dispersion}{
a numeric variable specifying the dispersion of the fit to be assumed when computing the standard errors.
}
  \item{terms}{
a character vector specifies which terms are to be returned. This is the case for \code{type = "terms"}.
}
  \item{na.action}{
function determining what should be done with missing values in \code{newdata}. The default is to predict \code{NA}.
}
  \item{\dots}{
other arguments that are passed into the function.
}
}
\details{
The specifications of the arugments are the same as the setting in \code{glm} function. The user may refer to the \code{\link[stats]{predict.glm}}.
}
\author{
Qihuang Zhang and Grace Y. Yi.
}

\seealso{
\code{\link[stats]{predict}},\code{\link[stats]{predict.glm}}
}

\examples{
data(ToyUni)
example <- augSIMEX(mainformula = Y ~ Xstar + Zstar + W, family = binomial(link = logit),
  mismodel = pi|qi ~ W, 
  meformula = Xstar ~ X + Z + W,
  data = ToyUni$Main,validationdata = ToyUni$Validation, subset = NULL,
  err.var = "Xstar", mis.var = "Zstar", err.true = "X", mis.true = "Z", 
  err.mat = NULL,
  lambda = NULL, M = 5, B = 2, nBoot = 2, extrapolation="quadratic")
predict(example)
}

\keyword{function}
