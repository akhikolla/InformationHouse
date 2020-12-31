\name{augSIMEX}
\alias{augSIMEX}
\title{
Analysis of Data with Mixed Measurement Error and Misclassification in Covariates
}
\description{
Implementation of the SIMEX algorithm for data with mixed measurement error and misclassification in covariates.
}
\usage{
augSIMEX(mainformula = formula(data), mismodel = pi | qi ~ 1,
                 meformula = NULL, family = gaussian, data,
                 validationdata, err.var, mis.var, err.true, mis.true,
                 err.mat = NULL, cppmethod = TRUE, repeated = FALSE,
                 repind = list(), subset, offset, weights, na.action,
                 scorefunction = NULL, lambda = NULL, M = 5, B = 20,
                 nBoot = 50, extrapolation = c("quadratic", "linear"), bound = 8,
                 initial = NULL, ...)
}
\arguments{
  \item{mainformula}{
an object of class ``formula'': an object of class ``formula'': a symbolic description of the model of the response variable, error-prone covariates, and other covariates.
}
  \item{mismodel}{
an object of class ``Formula''. A symbolic description in modeling the misclassification rates. See details for the specification of the model.
}
  \item{meformula}{
an object of ``Formula'' specifying the measurement error model for each error-prone covariate. The number of responses should equal the number of error-prone covariates. The default choice will be classic additive models. See details for the specification of the model.
}
  \item{family}{
an object of class ``family'' (same as in \code{\link[stats]{glm}} function in stats package).
}
  \item{data}{
a data frame or a matrix of main data. The variable in the main data includes the response, observed covariates matrix which is subject to measurement error, observed binary covariate vector that is prone to misclassification and may also contain the precisely measured covariates matrix. The default choice includes all covariates that mentioned in the \code{mainformula}.
}
  \item{validationdata}{
a data frame or matrix of validation data. The variable in the model includes an observed covariates matrix of that is subject to measurement error and their corresponding precisely measured covariates, an observed binary covariate vector that is prone to misclassification and their corresponding precisely measured covariates, a covariates matrix of precisely measured covariates. 
}
  \item{err.var}{
a vector of character specifying the name of covariates that are subject to measurement error.
}
  \item{mis.var}{
a string specifying the variable name of the binary variable subject to misclassification.
}
  \item{err.true}{
a vector of character specifying the names of the all the counterparts of \code{err.var} that are precisely measured covariates. 
}
  \item{mis.true}{
a string specifying the variable name of the precisely measured binary variable.
}
  \item{err.mat}{
a matrix indicating the variance-covariance matrix of generated Normal distribution in simulation step. 
}
  \item{cppmethod}{
a logical value indicating whether solving the score function via C++ functions. The function involves \code{Rcpp} package. The C++ based method is much faster and more computationally efficient than the R based method. See the discussion section.
}
 \item{repeated}{
a logic value indicating whether repeated measurements are involved.
  }
  \item{repind}{
a list of vectors of repeated measurement for error-prone covariates. See detail below.
  }
  \item{subset}{
a logical vector indicating which subjects should be included in the fitting.
}
  \item{offset}{
a numeric vector indicating the offset. The default choice is null.
}
  \item{weights}{
an optional numeric vector of the weights should be used in the fitting.
}
  \item{na.action}{
a function indicating what method should be applied to deal with missing data.
}
  
  \item{scorefunction}{
a function of score function. To allow for the generality, the users can specify their own score function. It should be a function of parameters, response, covariates, weights and offset. A matrix of score value for each individual and each parameter should be returned in the function. An example is shown in \code{\link{glmscore}}
}
  \item{lambda}{
if the \code{M} value is not specified, the user can specify a positive sequence of lambda directly. The first element is usually set to be 0.
}
  \item{M}{
the number of predetermined lambda vector if \code{lambda} is missing. The default is chosen to be 50.
}
  \item{B}{
the number of the dataset generated in the simulation step. The default value is set to be 200.
}
  \item{nBoot}{
the number of the iterations repeated in bootstrap. The default value is set to be 100.
}
  \item{extrapolation}{
specifies the regression model that involves in the extrapolation step. The options include ``linear'' and ``quadratic''. The default is set to be ``quadratic''.
}
  \item{bound}{ a value or vector specififying the bound for the absolute value of the regression coefficients. During the simulation, the parameters out of bound will be filter out before extrapolation step.
  }
  \item{initial}{ the initial value of the parameters.
  }
  \item{\dots}{
other arguments that pass into the function.
}
}
\details{
The misclassification models are set in "Formula" format, where the misclassification rates for both classes of the binary variable is set simultaneously. The left-hand side sets the responses of the misclassification model. The response should always set be to pi (indicating pimodel) and qi (indicating qimodel), separated by "|". On the right hand side of the formula sets the covariates of the misclassification model. If the covariates for both models are different, use "|" to separate them in the same order as the response. See example.

The measurement error models are also set in "Formula" format. The left-hand side sets the responses of the measurement error model, which should be consistent with the specification of \code{err.var}. Each response is separated by "|". On the right-hand side of the formula, the covariates are set. If the covariates for each response are different, use "|" to separate the specifications of covariates and the order should correspond to that of the responses of the measurement error model. See example.

In the case of repeated measurements, the users can pick one of the measurements into the formula of the main model. If more than one covariates have multiple replicates, the users should name the vector of repeated measurements in the list of \code{repind} by the corresponding representative measurement in the formula of the main model. See example.

The number of bootstrap repetitions should be at least 2. i.e., \code{nBoot}>2. Otherwise, an error might occur. 

The examples are mainly for illustration purpose. NA's are possible to be generated because the bootstrap simulation parameter is only set as 2. To obtain precise results, simulation parameters are supposed to be set on a larger scale, which also involves more time for computing. 
}
\value{
\item{coefficients}{the coefficient of the main model after correction.}
\item{vcov}{a adjusted variance-covariance matrix estimated by bootstrap.}
\item{fitted.values}{the fitted mean values, obtained by transforming the linear predictors by the inverse of the link function.}
\item{lambda}{a vector of values that involved in extrapolation step.}
\item{coefmatrix}{a matrix of coefficient before doing extrapolation. Each row correspond to a specified lambda. Each column represent a component of coefficient. This output if for the plotting purpose.}
\item{deviance}{minus twice the maximized log-likelihood.}
\item{null.deviance}{the deviance of null model which only includes intercept and offset term. It is comparable to the deviance. }
\item{df.residual}{the degree of freedom.}
\item{df.null}{the degree of freedom of null model.}
\item{aic}{Akaike's Information Criterion, minus twice the maximized log-likelihood plus twice the number of parameters.}
\item{...}{Other components are the arguments that have been used in the function.}


}
\references{
Yi G Y, Ma Y, Spiegelman D, et al. Functional and structural methods with mixed measurement error and misclassification in covariates[J]. Journal of the American Statistical Association, 2015, 110(510): 681-696.
}
\author{
Qihuang Zhang and Grace Y. Yi
}

\seealso{
\code{\link[stats]{glm}}
}
\examples{
### Example 1: Univariate Case
data(ToyUni)
example1<-augSIMEX(mainformula = Y ~ Xstar + Zstar + W, family = binomial(link = logit),
  mismodel = pi|qi ~ W, 
  meformula = Xstar ~ X + Z + W,
  data = ToyUni$Main,validationdata = ToyUni$Validation, subset = NULL,
  err.var = "Xstar", mis.var = "Zstar", err.true = "X", mis.true = "Z", 
  err.mat = NULL,
  lambda = NULL, M = 5, B = 2, nBoot = 2, extrapolation="quadratic")                   
summary(example1)

## Without adjustment
example1_naive <- glm(formula = Y ~ Xstar + Zstar + W,
family = binomial(link = logit),data = ToyUni$Main)
summary(example1_naive)

## using  accurate data
example1_true <- glm(Y~Xstar+Zstar+W, family = binomial(link=logit), 
   data=ToyUni$True)


### Example 2: Multivariate Case
data(ToyMult)
ErrorFormula<-Xstar.X1|Xstar.X2~-1+X.X1|-1+X.X2   ## measurement error model
example2<-augSIMEX(mainformula = Y~Xstar.X1+Xstar.X2+Zstar+W.W1+W.W2, 
  mismodel=pi|qi~X.X1+X.X2+W.W1+W.W2, family = binomial,
  data=ToyMult$Main,meformula=ErrorFormula,
  validationdata=ToyMult$Validation, subset=NULL,
  err.var=c("Xstar.X1","Xstar.X2"), mis.var="Zstar", err.true=c("X.X1","X.X2"), 
  mis.true="Z", err.mat = NULL,
  lambda=NULL, M=5, B=2, nBoot=2, extrapolation="quadratic")
summary(example2)


### Example 3
data(ToyRepeat)
example3<-augSIMEX(mainformula = Y~Xstar1+Zstar+W, family = binomial(link=logit),
  mismodel = pi|qi ~ W, meformula = Xstar ~ X + Z + W,
  data=ToyRepeat$Main,validationdata=ToyRepeat$Validation, 
  subset=NULL, err.var="Xstar1", mis.var="Zstar", err.true="X", mis.true="Z", 
  err.mat = NULL, repeated = TRUE,repind=list(Xstar1=c("Xstar1","Xstar2")),
  lambda=NULL, M=5, B=2, nBoot=2, extrapolation="quadratic")
summary(example3)
}
\keyword{funtion}
\keyword{core}
