\name{predict.fk_ppr}
\alias{predict.fk_ppr}
\title{Predict method for class fk_ppr}
\description{
  Standard prediction method for regression models, specific to outputs from the fk_ppr() function. See help(fk_ppr) for more details on the model.
}
\usage{
\method{predict}{fk_ppr}(object, Xtest = NULL, ...)
}
\arguments{
  \item{object}{an object of class fk_ppr, output from the function of the same name.}
  \item{Xtest}{(optional) matrix of test data on which predictions are to be made. If omitted then
  fitted values from training data are returned.}
  \item{...}{(optional) further arguments passed to or from other methods.}
}
\value{
  A vector of predictions for Xtest.
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

### estimate comparative estimated coefficients of determination

MSE1 = mean((predict(model1, X[-train_ids,])-y[-train_ids])^2)
MSE2 = mean((predict(model2, X[-train_ids,])-y[-train_ids])^2)
Var = mean((y[-train_ids]-mean(y[-train_ids]))^2)


1-MSE1/Var
1-MSE2/Var
}
\keyword{file}
