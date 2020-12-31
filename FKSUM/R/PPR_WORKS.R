# f_ppr evaluates the projection index for projection pursuit regression.
# Arguments:
# v = projection vector
# X = matrix of covariates
# y = vector of responses
# h = positive numeric bandwidth
# betas = vector of kernel coefficients
# loss = loss function to be used. Must be a function which takes two arguments: a vector of responses and a vector of fitted values and returns
#         a vector containing the losses/costs for pairs of response and fitted value
# dloss = derivative of loss function. Must be a function which takes the same arguments as loss, and returns the vector of partial derivatives
#         of the elements in loss(y, yhat) with respect to the elements in yhat. Not used in function evaluation, but in place so that arguments match with df_ppr
# type = character. One of 'loc-lin' for local linear regression and 'NW' for Nadaraya-Watson (local constant) regression
#
# returns a numeric, the output from loss(y, yhat)

f_ppr <- function(v, X, y, h, betas, loss = NULL, dloss = NULL, type = 'loc-lin'){
  # evaluate projected covariates and set-up
  p <- X %*% v / sqrt(sum(v^2))
  n <- length(p)
  o <- order(p)
  po <- p[o]
  yo <- y[o]

  # the fitted values differ depending on the type of regression estimator
  if(type == 'loc-lin'){ # local linear regression
    # compute all kernel sums needed to determine fitted values
    sK <- ksum(po, rep(1, n), po, h, betas, 1:n) - betas[1]
    sKx <- ksum(po, po, po, h, betas, 1:n) - po * betas[1]
    sKy <- ksum(po, yo, po, h, betas, 1:n) - yo * betas[1]
    sKx2 <- ksum(po, po^2, po, h, betas, 1:n) - po^2 * betas[1]
    sKxy <- ksum(po, po * yo, po, h, betas, 1:n) - po * yo * betas[1]

    # evaluate fitted values
    hy <- ((sKx2 * sKy - sKx * sKxy) + (sK * sKxy - sKx * sKy) * po) / (sK * sKx2 - sKx^2)
  }
  else{ # Nadaraya-Watson regression
    # compute all kernel sums needed to determine fitted values
    sK <- ksum(po, rep(1, n), po, h, betas, 1:n) - betas[1]
    # buffer values away from zero for more stable performance
    sK[sK < 1e-20] <- 1e-20
    sKy <- ksum(po, yo, po, h, betas, 1:n) - yo * betas[1]

    # evaluate fitted values
    hy <- sKy / sK
  }
  # compute the loss for the current projection
  if(is.null(loss)) sum((yo - hy)^2)
  else sum(loss(yo, hy))
}

# kLLreg evaluates the fitted values, almost as fk_regression(), except maintaining the order of the inputs.
# Operation of the function is exactly as the fitting of hy in f_ppr, but without leave-one-out sums.
# Arguments:
# x = vector of univariate covariates. In PPR context, the projected data X %*% v
# y = vector of response values
# h = positive numeric bandwidth value
# betas = vector of kernel coefficients
# type = type of regression to use, as in f_ppr, etc.

kLLreg <- function(x, y, h, betas, type){
  n <- length(x)
  o <- order(x)
  xo <- x[o]
  yo <- y[o]
  if(type == 'loc-lin'){
    sK <- ksum(xo, numeric(n) + 1, xo, h, betas)
    sKx <- ksum(xo, xo, xo, h, betas)
    sKy <- ksum(xo, yo, xo, h, betas)
    sKx2 <- ksum(xo, xo^2, xo, h, betas)
    sKxy <- ksum(xo, xo * yo, xo, h, betas)
    (((sKx2 * sKy - sKx * sKxy) + (sK * sKxy - sKx * sKy) * xo) / (sK * sKx2 - sKx^2))[rank(x)]
  }
  else{
    sK <- ksum(xo, numeric(n) + 1, xo, h, betas)
    sKy <- ksum(xo, yo, xo, h, betas)
    (sKy / sK)[rank(x)]
  }
}


# df_ppr evaluates the gradient of the projection index for projection pursuit regression.
# Arguments:
# v = projection vector
# X = matrix of covariates
# y = vector of responses
# h = positive numeric bandwidth
# betas = vector of kernel coefficients
# loss = loss function to be used. Must be a function which takes two arguments: a vector of responses and a vector of fitted values and returns
#         a vector containing the losses/costs for pairs of response and fitted value
# dloss = derivative of loss function. Must be a function which takes the same arguments as loss, and returns the vector of partial derivatives
#         of the elements in loss(y, yhat) with respect to the elements in yhat. Not used in function evaluation, but in place so that arguments match with df_ppr
# type = character. One of 'loc-lin' for local linear regression and 'NW' for Nadaraya-Watson (local constant) regression
#
# returns a vector, the partial derivatives of f_ppr(v, ...) with respect to the elements in v

df_ppr <- function(v, X, y, h, betas, loss = NULL, dloss = NULL, type = 'loc-lin'){
  # compute projected covariates and set-up
  nv <- sqrt(sum(v^2))
  p <- X %*% v / nv
  n <- length(p)
  o <- order(p)
  po <- p[o]
  yo <- y[o]

  # the fitted values, and hence their partial derivatives, differ depending on the type of regression estimator
  if(type == 'loc-lin'){ # local linear regression
    # compute all kernel (derivative) sums needed to evaluate the partial derivatives of f_ppr(v, X, ...) w.r.t.
    # the elements in X %*% v / norm(v)
    S1 <- kndksum(po, numeric(n) + 1, po, h, betas, 1:n) - rep(c(betas[1], 0), each = n)
    Sx <- kndksum(po, po, po, h, betas, 1:n) - cbind(po * betas[1], 0)
    Sy <- kndksum(po, yo, po, h, betas, 1:n) - cbind(yo * betas[1], 0)
    Sxx <- kndksum(po, po^2, po, h, betas, 1:n) - cbind(po^2 * betas[1], 0)
    Sxy <- kndksum(po, po * yo, po, h, betas, 1:n) - cbind(po * yo * betas[1], 0)
    N <- ((Sxx[, 1] * Sy[, 1] - Sx[, 1] * Sxy[, 1]) + (S1[, 1] * Sxy[, 1] - Sx[, 1] * Sy[, 1]) * po)
    D <- (S1[, 1] * Sxx[, 1] - Sx[, 1]^2)
    hy <- N / D
    if(is.null(loss)){
      dL <- 2 * (hy - yo)
    }
    else dL <- dloss(yo, hy)
    dLD <- dL / D
    Sy.p_Spy_Sp.hy <- (Sy[, 1] * po + Sxy[, 1] - 2 * Sx[, 1] * hy) * dLD
    S1.hy_Sy <- (S1[, 1] * hy - Sy[, 1]) * dLD
    Sp_S1.p <- (Sx[, 1] - S1[, 1] * po) * dLD
    Sp.p_Spp <- (Sx[, 1] * po - Sxx[, 1]) * dLD
    Spp.hy_Spy.p <- (Sxx[, 1] * hy - Sxy[, 1] * po) * dLD

    S_Sy.p_Spy_Sp.hy <- kndksum(po, Sy.p_Spy_Sp.hy, po, h, betas, 1:n) - betas[1] * cbind(Sy.p_Spy_Sp.hy, 0)
    S_S1.hy_Sy <- kndksum(po, S1.hy_Sy, po, h, betas, 1:n) - betas[1] * cbind(S1.hy_Sy, 0)
    S_Sp_S1.p <- kndksum(po, Sp_S1.p, po, h, betas, 1:n) - betas[1] * cbind(Sp_S1.p, 0)
    S_Sp.p_Spp <- dksum(po, Sp.p_Spp, po, h, betas, 1:n)
    S_Spp.hy_Spy.p <- dksum(po, Spp.hy_Spy.p, po, h, betas, 1:n)

    # compute the vector of partial derivatives of f_ppr(v, X, ...) w.r.t. the elements in X %*% v / norm(v)
    dp <- po / h * S_Sy.p_Spy_Sp.hy[, 2] + po^2 / h * S_S1.hy_Sy[, 2] + po * yo / h * S_Sp_S1.p[, 2] + yo / h * S_Sp.p_Spp
    dp <- dp + S_Spp.hy_Spy.p / h - 2 * po * S_S1.hy_Sy[, 1] - S_Sy.p_Spy_Sp.hy[, 1] - yo * S_Sp_S1.p[, 1]
    dp <- dp + dLD * (Sy[, 1] * (po / h * Sx[, 2] - 1 / h * Sxx[, 2] - Sx[, 1]) - 1 / h * Sy[, 2] * (Sxx[, 1] - po * Sx[, 1]) - 1 / h * Sxy[, 2] * (po * S1[, 1] - Sx[, 1]))
    dp <- dp + dLD * (Sxy[, 1] * (Sx[, 2] / h - po / h * S1[, 2] + S1[, 1]) + hy * (Sxx[, 1] / h * S1[, 2] + Sxx[, 2] / h * S1[, 1] - 2 / h * Sx[, 1] * Sx[, 2]))
  }
  else{ # Nadaraya-Watson (local constant) regression
    # compute all kernel (derivative) sums needed to evaluate the partial derivatives of f_ppr(v, X, ...) w.r.t.
    # the elements in X %*% v / norm(v)
    S1 <- kndksum(po, numeric(n) + 1, po, h, betas, 1:n) - rep(c(betas[1], 0), each = n)
    S1[S1[, 1] < 1e-20, 1] <- 1e-20
    Sy <- kndksum(po, yo, po, h, betas, 1:n) - cbind(yo * betas[1], 0)
    hy <- Sy[, 1] / S1[, 1]
    if(is.null(loss)){
      dL <- 2 * (hy - yo)
    }
    else dL <- dloss(yo, hy)

    # compute the vector of partial derivatives of f_ppr(v, X, ...) w.r.t. the elements in X %*% v / norm(v)
    dp <- (dksum(po, hy * dL / S1[, 1], po, h, betas, 1:n) - yo * dksum(po, dL / S1[, 1], po, h, betas, 1:n) +
             dL / S1[, 1] * (hy * S1[, 2] - Sy[, 2])) / h
  }

  # return the vector of partial derivatives of f_ppr(v, ...) w.r.t. the elements of v, based on a chain rule product
  c(c(dp) %*% (X[o, ] / nv - po %*% t(v) / nv^2))
}


# fk_ppr performs projection pursuit regression. The function uses a forward fitting procedure which adds components which minimise
# the loss/error based on the residuals remaining after the components already in the model have been accounted for.
# Arguments:
# X = matrix of covariates
# y = vector of responses
# nterms = number of components to add to the model
# hmult = positive numeric scaling to apply to the default bandwidth rule used
# loss = the loss function to be passed to f_ppr and df_ppr
# dloss = the function to evaluate the derivative of the loss function, to be passed to f_ppr and df_ppr
# initialisation = the method to be used for initialisation of projection pursuit. Must be either a function of only the
#       covariates and responses (function(X, y)), which returns a vector of length ncol(X), or one of 'lm' to initialise
#       with linear regression coefficients, or 'random' for random initialisation.
# type = character. One of 'loc-lin' for local linear regression and 'NW' for Nadaraya-Watson (local constant) regression
#
# returns a named list with class "fk_ppr" with fields
# $mu = mean of the responses
# $mu_X = vector of means of the covariates
# $X, $y, $betas, $type = the corresponding arguments passed to the function
# $hs = vector of bandwidths used in each of the components of the model
# $vs = matrix whose rows are the optimal projection vectors from each of the components
# $fitted = matrix whose rows are the fitted values along each of the optimal projections, based on the residuals remaining at each step


fk_ppr <- function(X, y, nterms = 1, hmult = 1, betas = NULL, loss = NULL, dloss = NULL, initialisation = 'lm', type = 'loc-lin'){
  # check inputs
  if(!is.matrix(X)) stop('X must be a numeric matrix')
  if(!is.numeric(y) || length(y)!=nrow(X)) stop('y should be a numeric vector of length nrow(X)')
  if(any(is.na(c(X, y)))) stop('X and y cannot contain missing values')
  if(!is.null(betas) && (!is.vector(betas) || !is.numeric(betas))) stop('betas must be a numeric vector')
  if(!is.numeric(hmult) || length(hmult)>1 || hmult<0) stop('hmult must be a positive numeric')
  if(!is.null(loss) && !is.function(loss)) stop('loss must be a function to compute the losses from the pairs in y and yhat')
  if(!is.null(dloss) && !is.function(dloss)) stop('dloss must be a function to compute the partial derivatives of the total loss w.r.t. yhat, from arguments y and yhat')
  if(!initialisation%in%c('lm', 'random') && !is.function(initialisation)) stop('initialisation must be one of "lm" and "random" or a function of X and y returning a vector of length ncol(X)')
  if(!type%in%c('loc-lin', 'NW')) stop('type must be one of "loc-lin" and "NW"')

  # set-up and compute necessary parameters
  n <- nrow(X)
  d <- ncol(X)
  mu <- mean(y)
  mu_X <- colMeans(X)

  # r will be used to represent the residuals remaining for fitting each component
  # the first component is fit to the centralised response values
  r <- y - mu

  # resids is used to store the residuals used to obtain each component

  resids <- matrix(0, nterms, n)

  # computation is slightly more stable for covariates of smaller relative magnitude
  X <- sweep(X, 2, mu_X, '-')

  # set-up objects to store parts of the solution
  hs <- numeric(nterms)
  vs <- matrix(0, nterms, d)
  fitted <- matrix(0, nterms, n)
  if(is.null(betas)) betas = c(0.25, 0.25)

  # add the components in the model iteratively
  for(tm in 1:nterms){
    # determine initial projection vector
    if(initialisation == 'lm') v0 <- solve(t(X) %*% X + .01 * diag(ncol(X))) %*% t(X) %*% r
    else if(initialisation == 'random') v0 <- rnorm(d)
    else if(is.function(initialisation)) v0 <- initialisation(X, r)
    else stop('argument "initialisation" must be a function of X and y or one of "lm" and "random".')

    # normalise the initial projection
    v0 <- v0 / sqrt(sum(v0^2))

    # compute the bandwidth to be used during projection pursuit. We use an over-smoothing bandwidth at this stage for
    # more stable performance.
    if(d < 100) h <- sqrt(eigen(cov(X))$values[1]) / n^.2 * hmult
    else h <- sqrt(eigs_sym(cov(X), 1)$values[1]) / n^.2 * hmult

    # optimise the projection vector using the base function optim(), and the L-BFGS algorithm
    vs[tm,] <- optim(v0, f_ppr, df_ppr, X, r, h, betas, loss, dloss, type, method = 'L-BFGS-B', control = list(maxit = 50))$par

    # normalise the projection vector
    vs[tm,] <- vs[tm,] / sqrt(sum(vs[tm,]^2))

    # determine the bandwidth for use in the final fitted values based on cross-validation
    p <- X %*% vs[tm,]
    cv_fun <- function(h) f_ppr(matrix(1, 1, 1), p, r, h, betas, loss, dloss, type)
    h <- optimise(cv_fun, c(.05 * h, h))$minimum

    # compute the fitted values for the current component
    fitted[tm,] <- kLLreg(X %*% vs[tm,], r, h, betas, type)

    # store the residuals used to obtain the projection
    resids[tm,] <- r

    # update the residuals to be used for the subsequent component
    r <- r - fitted[tm,]

    hs[tm] <- h
  }
  # return the list with all the parts of the model, and necessary parameters/arguments, and set its class
  structure(list(mu = mu, mu_X = mu_X, y = y, X = X, hs = hs, vs = vs, fitted = fitted, res = resids, betas = betas, type = type, nterms = nterms, call = match.call()), class = "fk_ppr")
}


# Methods for class fk_ppr:



# predict.fk_ppr is the S3 method for the class "fk_ppr", and determines the predictions from a model output from
# the function fk_ppr on test data
# Arguments:
# object = a list with class "fk_ppr"; the output from the function fk_ppr
# Xtest = a matrix of test data (covariates) whose response values need to be estimated/predicted. If omitted
#       then the function returns the fitted values from the training data, X, passed to the function fk_ppr
# ... = not currently utilised
#
# returns a vector of predictions, one for each of the rows in Xtest.

predict.fk_ppr <- function(object, Xtest = NULL, ...){
  # set-up and extract some of the elements of object for brevity

  if(is.null(Xtest)) Xtest <- object$X
  else Xtest <- sweep(Xtest, 2, object$mu_X, '-')
  betas <- object$betas
  n <- nrow(object$X)
  ntest <- nrow(Xtest)

  # yhat stores the predictions, which begin with the mean of the responses from the training data
  yhat <- object$mu

  # the fitted values from each component in the model are now used to incrementally adjust the predictions
  # by interpolation using the appropriate regression method
  for(tm in 1:length(object$hs)){
    # determine the kernel locations (projected training data on the current component projection) and bandwidth
    h <- object$hs[tm]
    p <- object$X %*% object$vs[tm,]
    o <- order(p)

    # evaluate the regression at the projected test data, ptest
    ptest <- Xtest %*% object$vs[tm,]
    otest <- order(ptest)

    # the evaluation of the regression for each component is as in the function fk_regression
    if(object$type == 'loc-lin'){
      sK <- ksum(p[o], numeric(n) + 1, ptest[otest], h, betas)
      sKx <- ksum(p[o], p[o], ptest[otest], h, betas)
      sKy <- ksum(p[o], object$res[tm, o], ptest[otest], h, betas)
      sKx2 <- ksum(p[o], p[o]^2, ptest[otest], h, betas)
      sKxy <- ksum(p[o], p[o] * object$res[tm, o], ptest[otest], h, betas)
      yhat <- yhat + (((sKx2 * sKy - sKx * sKxy) + (sK * sKxy - sKx * sKy) * ptest[otest]) / (sK * sKx2 - sKx^2))[rank(ptest)]
    }
    else{
      sK <- ksum(p[o], numeric(n) + 1, ptest[otest], h, betas)
      sKy <- ksum(p[o], object$res[tm,o], ptest[otest], h, betas)
      yhat <- yhat + (sKy / sK)[rank(ptest)]
    }
  }

  # return the predicted values
  yhat
}


plot.fk_ppr <- function(x, term = NULL, ...){
  if(is.null(term)){
    plot(colSums(x$fitted) + mean(x$y), x$y - colSums(x$fitted) - mean(x$y), main = 'Residuals vs. Fitted Values', xlab = 'Fitted Values', ylab = 'Residuals', ...)
  }
  else{
    plot(x$X%*%x$vs[term,], x$res[term,], main = paste('Term ', term, sep = ''), xlab = "x'v", ylab = "Fitted Value", ...)
    points(x$X%*%x$vs[term,], x$fitted[term,], col = 2)
  }
}

print.fk_ppr <- function(x, ...){
  cat('Call: \n \n')
  print(x$call)
  cat('\nGoodness of fit: \n')
  cat(paste(x$nterms, 'terms \n'))
  cat(paste('RSS = ', round(sum((x$y - colSums(x$fitted) - mean(x$y))^2), 4), '\n'))
}
