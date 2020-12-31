\name{pencoxfrail}
\alias{pencoxfrail}
\alias{pencoxfrail}
\docType{package}
\title{
Regularization in Cox Frailty Models.}
\description{A regularization approach for Cox Frailty Models by penalization methods is provided.
}
\details{
The \code{pencoxfrail} algorithm is designed to investigate
the effect structure in the Cox frailty model, which is a
widely used model that accounts for heterogeneity in survival data.
Since in survival models one has to account for possible variation of
the effect strength over time the selection of the relevant features distinguishes between the folllowing cases:
covariates can have time-varying effects, can have time-constant effects
or be irrelevant. For this purpose, the following specific penality is applied on the vectors of B-spline coefficients \eqn{\alpha_k}, assuming \eqn{k=1,...,r} different, potentially time-varying effects, each expanded in \eqn{M} B-spline basis functions:
\deqn{\xi \cdot J(\zeta,\alpha) = \xi\left( \zeta\sum_{k=1}^{r} \psi\, w_{\Delta,k} ||\Delta_{M}\alpha_{k}||_{2} + (1-\zeta)\sum_{k=1}^{r} \phi\, v_{k} ||\alpha_{k}||_{2}\right)}{\xi*J(\zeta,\alpha) = \xi * \{ \zeta * \sum_k \psi * w_k * ||\Delta_M*\alpha_k||_2 + (1-\zeta) * \sum_k \phi* v_k * ||\alpha_k||_2 \}.}
This penalty is able to distinguish between these
types of effects to obtain a sparse representation that includes the relevant effects in a proper form.

The penalty is depending on two tuning parameters, \eqn{\xi} and \eqn{\zeta}, which have to be determined by a suitable technique, e.g. by 
(2-dimensional) K-fold cross validation.

The first term of the penalty controls the smoothness of the time-varying covariate effects, whereby for
values of \eqn{\xi} and \eqn{\zeta} large enough, all differences \eqn{(\alpha_{k,l} - \alpha_{k,l-1}),  l=2,... ,M}{(\alpha_k,l - \alpha_k,l-1),  l=2,... ,M}, are removed from the model, resulting in constant covariate effects.
As the B-splines of each variable with varying coefficients sum up to one,
a constant effect is obtained if all spline coefficients are set equal.
Hence, the first penalty term does not affect the spline's global level.
The second term penalizes all spline coefficients belonging to a single time-varying effect in the way
of a group LASSO and, hence, controls the selection of covariates.


\tabular{ll}{
Package: \tab pencoxfrail\cr
Type: \tab Package\cr
Version: \tab 1.0.1\cr
Date: \tab 2016-05-06\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
for loading a dataset type data(nameofdataset)
}

\usage{
pencoxfrail(fix=formula, rnd=formula, vary.coef=formula, data, xi, 
              adaptive.weights = NULL, control = list())
}     
\arguments{
  \item{fix}{a two-sided linear formula object describing the unpenalized
    fixed (time-constant) effects part of the model, with the response on the left of a
    \code{~} operator and the terms, separated by \code{+} operators, on
    the right. The response must be a survival object as returned by the \code{\link[survival]{Surv}} function.}  
  \item{rnd}{a two-sided linear formula object describing the
    random-effects part of the model, with the grouping factor on the left of a
    \code{~} operator and the random terms, separated by \code{+} operators, on
    the right.}
  \item{vary.coef}{a one-sided linear formula object describing the
    time-varying effects part of the model, with the time-varying terms, separated by \code{+} operators,
    on the right side of a \code{~} operator.}
  \item{data}{the data frame containing the variables named in the three preceding
    \code{formula} arguments.}
  \item{xi}{the overall penalty parameter that controls the strenght of both penalty terms in \eqn{\xi\cdot J(\zeta,\alpha)}{\xi*J(\zeta,\alpha)} and, hence, controls the overall amount of smoothness (up to constant effects) and variable selection for a given proportion \eqn{\zeta}.
  The optimal penalty parameter is a tuning parameter of the procedure that has to be determined, 
 e.g. by K-fold cross validation. (See details or the quick demo for an example.)} 
  \item{adaptive.weights}{a two-column matrix of adaptive weights passed to the procedure; the first column contains the weights \eqn{w_{\Delta,k}}{w_k}, the second column the weights \eqn{v_{k}}{v_k} from \eqn{\xi\cdot J(\zeta,\alpha)}{\xi*J(\zeta,\alpha)}. If no adaptive weights are specified all weights are set to one. The recommended strategy is to first fit an unpenalized model (i.e. \eqn{\xi=0}) and then use the obtained adaptive weights (see value section) when fitting the model for all other combinations of \eqn{\xi} and \eqn{\zeta}.} 
  \item{control}{a list of control values for the estimation algorithm to replace the default values returned by the function \code{\link{pencoxfrailControl}}. Defaults to an empty list.}
}


\value{Generic functions such as \code{print}, \code{predict}, \code{plot} and \code{summary} have methods to show the results of the fit.

The \code{predict} function uses also estimates of random effects for prediction, if possible (i.e. for known subjects of the grouping factor). 
Either the survival stepfunction or the baseline hazard (not cumulative!) can be calculated by specifying one of two possible methods: \code{method=c("hazard","survival")}. By default, for each new subject in \code{new.data} an individual stepfunction is calculated on a pre-specified time grid, also accounting for covariate changes over time. Alternatively, for \code{new.data} a single vector of a specific (time-constant) covariate combination can be specified.

Usage:  \code{
predict(pencoxfrail.obj,new.data,time.grid,method=c("hazard","survival"))
}     


The \code{plot} function plots all time-varying effects, including the baseline hazard. 

   \item{call}{a list containing an image of the \code{pencoxfrail} call that produced the object.}  
     \item{baseline}{a vector containing the estimated B-spline coefficients of the baseline hazard.
     If the covariates corresponding to the time-varying effects are centered (and standardized, see \code{\link{pencoxfrailControl}}), the coefficients are transformed back to the original scale.} 
          \item{time.vary}{a vector containing the estimated B-spline coefficients of all time-varying effects.
     If the covariates corresponding to the time-varying effects are standardized (see \code{\link{pencoxfrailControl}}) 
     the coefficients are transformed back to the original scale.} 
  \item{coefficients}{a vector containing the estimated fixed effects.}
  \item{ranef}{a vector containing the estimated random effects.}
  \item{Q}{a scalar or matrix containing the estimates of the random effects standard deviation or variance-covariance parameters, respectively.}
  \item{Delta}{a matrix containing the estimates of fixed and random effects (columns) for each iteration (rows) of the main algorithm (i.e. before the final re-estimation step is performed, see details).}
  \item{Q_long}{a list containing the estimates of the random effects variance-covariance parameters for each iteration of the main algorithm.}
  \item{iter}{number of iterations until the main algorithm has converged.}
  \item{adaptive.weights}{If \eqn{\xi=0}, a two-column matrix of adaptive weights is calculated; the first column contains the weights \eqn{w_{\Delta,k}}{w_k}, the second column the weights \eqn{v_{k}}{v_k} from \eqn{\xi\cdot J(\zeta,\alpha)}{\xi*J(\zeta,\alpha)}. If \eqn{\xi>0},
  the adaptive weights that have been used in the function's argument are displayed.}
  \item{knots}{vector of knots used in the B-spline representation.}
  \item{Phi.big}{large B-spline design matrix corresponding to the baseline hazard and all time-varying effects. For the time-varying effects, the B-spline functions (as a function of time) have already been multiplied with their associated covariates.}
  \item{time.grid}{the time grid used in when approximating the (Riemann) integral involved in the model's full likelihood.}
  \item{m}{number of metric covariates with time-varying effects.}
  \item{m2}{number of categorical covariates with time-varying effects.}
}



\author{
Andreas Groll  \email{groll@math.lmu.de}
}

\references{
Groll, A., T. Hastie and G. Tutz (2016). 
Regularization in Cox Frailty Models. Ludwig-Maximilians-University. \emph{Technical Report} 191.
}


\seealso{
\code{\link{pencoxfrailControl},\link[survival]{Surv},\link[survival]{pbc}}
}
\examples{
\dontrun{
data(lung)

# remove NAs
lung <- lung[!is.na(lung$inst),]

# transform inst into factor variable
lung$inst <- as.factor(lung$inst)

# Random institutional effect
fix.form <- as.formula("Surv(time, status) ~ 1")
vary.coef <- as.formula("~ age")

pen.obj <- pencoxfrail(fix=fix.form,vary.coef=vary.coef, rnd = list(inst=~1), 
              data=lung, xi=10,control=list(print.iter=TRUE))

# show fit
plot(pen.obj)

# predict survival curve of new subject, institution 1 and up to time 500
pred.obj <- predict(pen.obj,newdata=data.frame(inst=1,time=NA,status=NA,age=26),
              time.grid=seq(0,500,by=1))

# plot predicted hazard function
plot(pred.obj$time.grid,pred.obj$haz,type="l",xlab="time",ylab="hazard")

# plot predicted survival function
plot(pred.obj$time.grid,pred.obj$survival,type="l",xlab="time",ylab="survival")

# see also demo("pencoxfrail-pbc")
}}
\keyword{
Lasso, Shrinkage, Variable selection, Generalized linear mixed model
}


