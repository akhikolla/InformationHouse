\name{flam-package}
\alias{flam-package}
\docType{package}
\title{
Fit the Fused Lasso Additive Model
}
\description{
This package is called flam for "fused lasso additive model", which is proposed in Petersen, A., Witten, D., and Simon, N. (2014). Fused Lasso Additive Model. arXiv preprint arXiv:1409.5391. The fused lasso additive model provides an approach to fit an additive model in which each component is estimated to be piecewise constant with a small number of adaptively-chosen knots.

The main functions are: (1) \code{\link{flam}} and (2) \code{\link{flamCV}}. The first function \code{flam} fits the fused lasso additive model for a range of tuning parameters and provides the fits for all of these tuning parameters. The second function \code{flamCV} considers a range of tuning parameters and provides the fits, but also returns the optimal tuning parameters, as chosen using K-fold cross-validation.

}
\details{
\tabular{ll}{
Package: \tab flam\cr
Type: \tab Package\cr
Version: \tab 3.0\cr
Date: \tab 2015-07-26\cr
License: \tab GPL (>= 2)\cr
}
The package includes the following functions:
\code{\link{flam}}, \code{\link{flamCV}}, \code{\link{plot.flam}}, \code{\link{plot.flamCV}}, \code{\link{plot.flamSparsity}}, \code{\link{predict.flam}}, \code{\link{summary.flam}}, \code{\link{summary.flamCV}}, \code{\link{flamDOF}}, and \code{\link{sim.data}}.

}
\author{
Ashley Petersen

Maintainer: Ashley Petersen <ajpete@uw.edu>
}
\references{
Petersen, A., Witten, D., and Simon, N. (2014). Fused Lasso Additive Model. arXiv preprint arXiv:1409.5391.
}
\keyword{ package }

\examples{
#general example illustrating all functions
#see specific function help pages for details of using each function

#generate data
set.seed(1)
data <- sim.data(n = 50, scenario = 1, zerof = 10, noise = 1)

#fit model for a range of tuning parameters, lambda and alpha
#lambda sequence is chosen automatically if not specified
flam.out <- flam(x = data$x, y = data$y, alpha.seq = c(0.8, 0.9, 1))
#or fit model and select lambda using 2-fold cross-validation
#note: use larger 'n.fold' (e.g., 10) in practice
flamCV.out <- flamCV(x = data$x, y = data$y, alpha = 1, n.fold = 2)

#summarize all of the fits (in particular, the percent sparsity achieved)
#note: percent sparsity is the percentage of features estimated
#to have no relationship with outcome
summary(flam.out)
#the percent sparsity over the range of tuning parameters can also 
#be displayed in a plot
plot(summary(flam.out))
#or just summarize a single fit 
#we'll examine the fit with an index of 25. that is, lambda and alpha of
flam.out$all.lambda[25]; flam.out$all.alpha[25]
summary(flam.out, index = 25)
#lastly, we can summarize the fit chosen using cross-validation
summary(flamCV.out$flam.out, index = flamCV.out$index.cv)
#the lambda chosen by cross-validation is also available using
flamCV.out$lambda.cv

#plot the estimated relationships between predictors and outcome
#do this for a specific fit
plot(flam.out, index = 25)
#or for the fit chosen using cross-validation
plot(flamCV.out$flam.out, index = flamCV.out$index.cv)
#by default, up to 10 non-sparse features with the largest L2 norms will 
#be plotted, see '?plot.flam' for other optional function arguments

#this data is simulated, so with a little more work, we can compare the 
#true generating functions to the estimated function fits
#we do this for the truly non-zero functions (i.e., the first four predictors)
#generate data from same model but larger n, just used to plot true functions
temp.data <- sim.data(n = 500, scenario = 1, zerof = 10, noise = 1)
col.vec = c("dodgerblue1","orange","seagreen1","hotpink")
theta.hat = flamCV.out$flam.out$theta.hat.list[[flamCV.out$index.cv]]
par(mfrow=c(2,2))
for (i in 1:4) {
	rgb.num = col2rgb(col.vec[i])
	col=rgb(rgb.num[1], rgb.num[2], rgb.num[3], 100, max=256)
	plot(1,type="n",xlim=c(-2.5,2.5),ylim=c(-2,2),xlab=paste("x",i,sep=""),
		ylab=paste("f",i,"(x",i,")",sep=""),main="")
	points(sort(temp.data$x[,i]), temp.data$theta[order(temp.data$x[,i]),i],type="l",lwd=3)
	points(sort(data$x[,i]), theta.hat[order(data$x[,i]),i],col=col,type="l",lwd=3)
}

#we can make predictions for a covariate matrix with new observations
#choose the alpha and lambda of interest
alpha <- flamCV.out$alpha; lambda <- flamCV.out$lambda.cv
#new.x with 20 observations and the same number of features as flam.out$x
new.data <- sim.data(n = 20, scenario = 1, zerof = 10, noise = 1)
new.x <- new.data$x
#these will give the same predictions:
yhat1 <- predict(flam.out, new.x = new.x, lambda = lambda, alpha = alpha)
yhat2 <- predict(flamCV.out$flam.out, new.x = new.x, lambda = lambda, alpha = alpha)

#we can summarize the cross-validation function call
summary(flamCV.out)
#and also plot the cross-validation error
plot(flamCV.out)
#or calculate degrees of freedom for the model chosen using cross-validation
flamDOF(object = flamCV.out$flam.out, index = flamCV.out$index.cv)
#or for any fit of a 'flam' object
flamDOF(object = flam.out, index = 25)
}
