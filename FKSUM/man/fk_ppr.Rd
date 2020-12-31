\name{fk_ppr}
\alias{fk_ppr}
\title{Projection pursuit regression with local linear kernel smoother}
\description{
  Generates a projection pursuit regression model for covariates X and response y
}
\usage{
fk_ppr(X, y, nterms=1, hmult=1, betas=NULL, loss=NULL,
    dloss=NULL, initialisation="lm", type = "loc-lin")
}
\arguments{
  \item{X}{a numeric matrix (num_data x num_dimensions) of covariates.}
  \item{y}{a numeric vector of responses.}
  \item{nterms}{The number of terms in the additive regression model. The default is a single term.}
  \item{betas}{The coefficients in the expression of the kernel. See Hofmeyr (2019) for details. The
  default is the smooth order one kernel described in the paper.}
  \item{hmult}{The bandwidth in the kernel smoother is set to eigen(cov(X))values[1]^.5/nrow(X)^.2*hmult during projection pursuit. The final value is based on leave-one-out cross validation on the optimal projection.}
  \item{loss}{The (additive) loss function to be used. Allows for generalised responses by specifying an appropriate likelihood/deviance function for the loss. Note: loss is to be minimised, so deviance or negative log-likelihood would be appropriate. If specified then must be a function of y and hy (the fitted values, yhat), returning a vector of "errors"" which are added as the total loss. The default is the squared error.}
  \item{dloss}{The derivative of the loss function. Also takes arguments y and hy, and returns the vector of partial derivatives of the loss w.r.t. hy.}
  \item{initialisation}{Method use to initialise the projection vectors. Must be one of "lm" and "random", or a function taking only arguments X and y, and returning a vector of length ncol(X). The default is "lm", which is a simple linear model with a small ridge to ensure a solution. "random" uses random initialisation.}
  \item{type}{The type of regression estimator. Must be one of "loc-lin" for local linear regression, or "NW" for the Nadaray
  Watson, or local constant regression. The default is local linear regression.}
}
\value{
  A named list with class fk_ppr containing
  \item{$mu}{The estimated (global) mean of the response.}
  \item{$mu_X}{The vector of estimated means of the covariates.}
  \item{$y}{The responses given as argument to the function.}
  \item{$X}{The covariates given as argument to the function.}
  \item{$hs}{A vector of bandwidths used for each term in the model.}
  \item{$vs}{A matrix whose rows are the projection vectors.}
  \item{$fitted}{The fitted values on each projection. Note that these are based on the residuals used to fit that component of the model, and not the original y values. $fitted is used for prediction on test data.}
  \item{$beta}{The beta coefficients in the kernel formulation.}
}
\references{
  Friedman, J., and Stuetzle, W. (1981) "Projection pursuit regression." \emph{Journal of the American statistical Association} \bold{76}.376.

Hofmeyr, D.P. (2019) "Fast exact evaluation of univariate kernel sums", \emph{IEEE Transactions on Pattern Analysis and Machine Intelligence}, in press.
}
\examples{
op <- par(no.readonly = TRUE)

set.seed(2)

### Generate a set of data

X = matrix(rnorm(10000), ncol = 10)

### Generate some "true" projection vectors

beta1 = (runif(10)>.5)*rnorm(10)
beta2 = (runif(10)>.5)*rnorm(10)

### Generate responses, dependent on X\%*\%beta1 and X\%*\%beta2

y = 1 + X\%*\%beta1 + ((X\%*\%beta2)>2)*(X\%*\%beta2-2)*10
y = y + (X\%*\%beta1)*(X\%*\%beta2)/5 + rnorm(1000)

### Fit a PPR model with 2 terms on a sample of the data

train_ids = sample(1:1000, 500)

model = fk_ppr(X[train_ids,], y[train_ids], nterms = 2)

### Predict on left out data, and compute
### estimated coefficient of determination

yhat = predict(model, X[-train_ids,])

MSE = mean((yhat-y[-train_ids])^2)
Var = mean((y[-train_ids]-mean(y[-train_ids]))^2)

1-MSE/Var


#################################################

### Add some "outliers" in the training data and fit
### the model again, as well as one with an absolute loss

y[train_ids] = y[train_ids] + (runif(500)<.05)*(rnorm(500)*100)

model1 <- fk_ppr(X[train_ids,], y[train_ids], nterms = 2)

model2 <- fk_ppr(X[train_ids,], y[train_ids], nterms = 2,
    loss = function(y, hy) abs(y-hy),
    dloss = function(y, hy) sign(hy-y))

### Plot the resulting components in the model on the test data

par(mar = c(2, 2, 2, 2))
par(mfrow = c(2, 2))

plot(X[-train_ids,]\%*\%model1$vs[1,], y[-train_ids])
plot(X[-train_ids,]\%*\%model1$vs[2,], y[-train_ids])
plot(X[-train_ids,]\%*\%model2$vs[1,], y[-train_ids])
plot(X[-train_ids,]\%*\%model2$vs[2,], y[-train_ids])

par(op)

### Estimate comparative estimated coefficients of determination

MSE1 = mean((predict(model1, X[-train_ids,])-y[-train_ids])^2)
MSE2 = mean((predict(model2, X[-train_ids,])-y[-train_ids])^2)
Var = mean((y[-train_ids]-mean(y[-train_ids]))^2)


1-MSE1/Var
1-MSE2/Var
}
\keyword{file}
