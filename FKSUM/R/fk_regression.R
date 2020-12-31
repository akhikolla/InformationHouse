# fk_regression computes kernel based estimate of regression function f(x) = E[Y|X=x]
# Arguments:
# x = numeric (univariate) covariate. The "independent variable" in some terminology
# y = numeric response variable
# h = bandwidth for use in the kernel estimate. Can be positive numeric or one of 'amise' or 'cv'. default is 'amise'
#         'amise' computes a rough estimate of the asymptotic mean integrated squared error minimiser.
#         'cv' computes h based on leave-one-out squared error cross validation.
# beta = vector of kernel coefficients. default is c(.25, .25)
# from, to = min and max of the points at which to evaluate the regression function if grid evaluation is desired.
#         default is min(x), max(x)
# ngrid = number of grid points at which to evaluate regression function (if grid evaluation is desired). default is 1000.
#         if both ngrid and nbin set to NULL then exact evaluation at the sample points, x, themselves is returned.
# nbin = number of bins to use if binned approximation is desired. default is exact evaluation on a grid of 1000 points
# type =  one of 'loc-lin' for local linear estimate and 'NW' for Nadaraya-Watson (local constant) estimator. default is 'loc-lin'
#
# returns a list with components
# $x = points at which regression function is estimated. If the original samples, then the output will be ordered
#         according to the order statistics of x.
# $y = value of estimated regression function
# $h = bandwidth value used in estimate


fk_regression <- function(x, y, h = 'amise', beta = NULL, from = NULL, to = NULL, ngrid = 1000, nbin = NULL, type = 'loc-lin'){
  # check inputs
  if(!is.numeric(x) || !is.numeric(y) || length(x)!=length(y)) stop('x and y must be numeric vectors of the same length')
  if(any(is.na(c(x, y)))) stop('x and y cannot contain missing values')
  if(!is.null(beta) && (!is.vector(beta) || !is.numeric(beta))) stop('beta must be a numeric vector')

  n <- length(x)

  # computation is slightly more stable for smaller magnitude kernel locations, so centralise x at zero
  mn <- mean(x)
  x <- x - mn
  if(!is.null(from)) from <- from - mn
  if(!is.null(to)) to <- to - mn

  # if exact evaluation performed then fast kernel computations require sorted data
  if(is.null(nbin) & is.unsorted(x)){
    o <- order(x)
    xo <-  x[o]
    yo <- y[o]
  }
  else{
    xo <- x
    yo <- y
  }
  if(is.null(beta)) beta <- c(.25, .25)

  # compute bandwidth
  if(h=='amise'){

    # AMISE bandwidth requires estimate of the curvature of the regression function.
    # we use an over-smoothing pilot bandwidth for this estimation
    h0 <- sd(x) / n^.2

    if(!is.null(nbin)){ # binning approximation used

      # estimate first derivative of the regressiong function, df
      # we use the NW estimator for this purpose for simplicity.
      # to obtain an estimate of the curvature we simply use the
      # first differences from the first derivative
      xs <- seq(min(x), max(x), length = nbin)
      wts1 <- sm_bin_wts(x, rep(1, n), nbin, xs[1], xs[nbin])
      wtsy <- sm_bin_wts(x, y, nbin, xs[1], xs[nbin])
      sK <- ksum(xs, wts1, xs, h0, beta, 1:nbin)
      sKy <- ksum(xs, wtsy, xs, h0, beta, 1:nbin)
      dsK <- dksum(xs, wts1, xs, h0, beta, 1:nbin)
      dsKy <- dksum(xs, wtsy, xs, h0, beta, 1:nbin)
      df <- (dsK * sKy - sK * dsKy) / sK^2 / h0

      # since binned approximation is used, we cannot directly estimate residual variance function
      # we use a very small bandwidth and compute the squared differences between the over- and
      # -undersmoothed functions to estimate this variance
      sK2 <- ksum(xs, wts1, xs, h0 / 1000, beta, 1:nbin)
      sKy2 <- ksum(xs, wtsy, xs, h0 / 1000, beta, 1:nbin)
      sig2 <- mean((sKy2 / sK2 - sKy / sK)^2)

      # compute estimate of AMISE bandwidth
      h <- (sig2 * roughness_K(beta) / var_K(beta)^2 / sum(sK[-1] / (n - 1) / h0 * diff(df)^2 / (xs[2] - xs[1])) / n / 2)^.2
    }
    else{ # exact evaluation of the estimate. The same steps are used as in the binned
          # approximation, except the estimate of the residual variance function is based
          # directly on (y-f(x))^2 values.
      sK <- ksum(xo, rep(1, n), xo, h0, beta, 1:n)
      sKy <- ksum(xo, yo, xo, h0, beta, 1:n)
      dsK <- dksum(xo, rep(1,n), xo, h0, beta, 1:n)
      dsKy <- dksum(xo, yo, xo, h0, beta, 1:n)
      df <- (dsK * sKy - sK * dsKy) / sK^2 / h0
      sig2 <- mean((yo - sKy / sK)^2)
      h <- (sig2 * roughness_K(beta) / var_K(beta)^2 / mean((diff(df) / diff(xo))^2) / n / 2)^.2
    }
  }
  else if(h=='cv'){

    # cross validation estimate is straightforward, except we don't estimate this for the binned
    # approximation as if the data set is very large and insufficient bins are used the performance
    # is poor.
    if(!is.null(nbin)){
      message('Cross validation bandwidth estimation not implemented for binned estimator.
          Switching to AMISE optimal estimate \n')
      return(fk_regression(x, y, h = 'amise', beta, from, to, ngrid, nbin, type))
    }
    else{

      # if exact evaluation is performed then cross validation is performed
      # by computing the leave-one-out estimates of the regression function
      # That is, argument loo is set to one in the regression computation.
      if(type=='loc-lin'){
        loo_cv <- function(h){
          sum((fk_loc_lin(xo, yo, h, beta, loo = 1) - yo)^2)
        }
      }
      else{
        loo_cv <- function(h){
          sum((fk_NW(xo, yo, h, beta, loo = 1) - yo)^2)
        }
      }

      # set h as the minimiser of the squared error from the leave-one-out
      # estimates of the regression function at the sample points themselves
      h <- optimise(loo_cv, sd(x) / n^.2 * c(.05, 2))$minimum
    }
  }
  else if(!is.numeric(h) || h < 0) stop('"h" must be either positive numeric or one of "cv" and "amise"')

  # In the remainder the differences in the output are in the evaluation points (either a grid of values in [from, to] or
  # if either ngrid or nbin is not NULL, or the ordered sample xo) and the arguments passed to the functions
  # fk_loc_lin or fk_NW, as these perform the actual computation of the regression estimates.

  if(!is.null(nbin)){ # binned approximation, hence x_eval is a grid and regression function estimated with nbin = nbin
    if(is.null(from)) from <- min(x)
    if(is.null(to)) to <- max(x)
    if(type=='loc-lin'){
      structure(list(x_eval = seq(from, to, length = nbin) + mn, y_fitted = fk_loc_lin(xo, yo, h, beta, nbin = nbin, from = from, to = to), x = x + mn, y = y, h = h, type = type, betas = beta, call = match.call()), class = 'fk_regression')
    }
    else if(type=='NW'){
      structure(list(x_eval = seq(from, to, length = nbin) + mn, y_fitted = fk_NW(xo, yo, h, beta, nbin = nbin, from = from, to = to), x = x + mn, y = y, h = h, type = type, betas = beta, call = match.call()), class = 'fk_regression')
    }
    else stop('type must be one of "loc-lin" and "NW"')
  }
  else if(!is.null(ngrid)){ # exact grid evaluation, hence x_eval is a grid and regression function estimated with nbin = NULL (default)
    if(is.null(from)) from <- xo[1]
    if(is.null(to)) to <- xo[n]
    if(type=='loc-lin'){
      structure(list(x_eval = seq(from, to, length = ngrid) + mn, y_fitted = fk_loc_lin(xo, yo, h, beta, ngrid = ngrid, from = from, to = to), x = x + mn, y = y, h = h, type = type, betas = beta, call = match.call()), class = 'fk_regression')
    }
    else if(type=='NW'){
      structure(list(x_eval = seq(from, to, length = ngrid) + mn, y_fitted = fk_NW(xo, yo, h, beta, ngrid = ngrid, from = from, to = to), x = x + mn, y = y, h = h, type = type, betas = beta, call = match.call()), class = 'fk_regression')
    }
    else stop('type must be one of "loc-lin" and "NW"')
  }
  else{ # exact evaluation at sample, hence both nbin and ngrid left as default (NULL) in regression computation
    if(type=='loc-lin'){
      structure(list(x_eval = xo + mn, y_fitted = fk_loc_lin(xo, yo, h, beta), x = x + mn, y = y, h = h, type = type, betas = beta, call = match.call()), class = 'fk_regression')
    }
    else if(type=="NW"){
      structure(list(x_eval = xo + mn, y_fitted = fk_NW(xo, yo, h, beta), x = x + mn, y = y, h = h, type = type, betas = beta, call = match.call()), class = 'fk_regression')
    }
    else stop('type must be one of "loc-lin" and "NW"')
  }
}

# fk_loc_lin computes the local linear regression estimate from sample (x, y)
# Arguments:
# x = (univariate) covariate values. If not using binned approximation then must be sorted non-decreasing
# y = response values. If not using binned approximation then must be ordered according to non-decreasing x-values
# h = positive numeric bandwidth value for kernel estimate
# beta = vector of kernel coefficients
# nbin = number of bins if binning approximation used. Default is exact evaluation
# ngrid = number of grid points if evaluation on a grid. Default is exact evaluation at sample points, x
# from, to = min and max point of evaluation if either grid evalution or binned approximation to be used
# loo = integer: 1 for leave-one-out estimate and 0 for standard estimate. loo = 1 not sensible for
#         grid estimates and not recommended for binned approximation
#
# returns vector of estimated regression function values at evaluation points

fk_loc_lin <- function(x, y, h, beta, nbin = NULL, ngrid = NULL, from = NULL, to = NULL, loo = 0){
  n <- length(x)
  if(!is.null(nbin)){ # binned approximation

    # compute bin locations
    if(is.null(from)) from <- x[1]
    if(is.null(to)) to <- x[n]
    xs <- seq(from, to, length = nbin)

    # compute smoothed bin weights for all terms required for regression estimate
    wts1 <- sm_bin_wts(x, rep(1, n), nbin, xs[1], xs[nbin])
    wtsx <- sm_bin_wts(x, x, nbin, xs[1], xs[nbin])
    wtsy <- sm_bin_wts(x, y, nbin, xs[1], xs[nbin])
    wtsx2 <- sm_bin_wts(x, x^2, nbin, xs[1], xs[nbin])
    wtsxy <- sm_bin_wts(x, x * y, nbin, xs[1], xs[nbin])

    # compute kernel weighted sums of all terms above
    sK <- ksum(xs, wts1, xs, h, beta, 1:nbin) - loo * beta[1]
    sKx <- ksum(xs, wtsx, xs, h, beta, 1:nbin) - loo * beta[1] * wtsx / wts1
    sKy <- ksum(xs, wtsy, xs, h, beta, 1:nbin) - loo * beta[1] * wtsy / wts1
    sKx2 <- ksum(xs, wtsx2, xs, h, beta, 1:nbin) - loo * beta[1] * wtsx2 / wts1
    sKxy <- ksum(xs, wtsxy, xs, h, beta, 1:nbin) - loo * beta[1] * wtsxy / wts1

    # return local-linear regression estimates
    (((sKx2 * sKy - sKx * sKxy) + (sK * sKxy - sKx * sKy) * xs) / (sK * sKx2 - sKx^2))
  }
  else if(!is.null(ngrid)){ # exact evaluation on a grid

    # compute evaluation points
    if(is.null(from)) from <- x[1]
    if(is.null(to)) to <- x[n]
    xs <- seq(from, to, length = ngrid)

    # compute kernel weighted sums of all terms required in local-linear estimate
    sK <- ksum(x, rep(1, n), xs, h, beta)
    sKx <- ksum(x, x, xs, h, beta)
    sKy <- ksum(x, y, xs, h, beta)
    sKx2 <- ksum(x, x^2, xs, h, beta)
    sKxy <- ksum(x, x * y, xs, h, beta)

    # return local-linear estimate of regression function
    (((sKx2 * sKy - sKx * sKxy) + (sK * sKxy - sKx * sKy) * xs) / (sK * sKx2 - sKx^2))
  }
  else{ # exact evaluation at sample points

    # compute kernel weighted sums of all terms required for regression estimate.
    # If loo = 1 then K(0)*(terms being summed) is subtracted from the weighted
    # sums to obtain leave-one-out estimates
    sK <- ksum(x, rep(1, n), x, h, beta, 1:n) - loo * beta[1]
    sKx <- ksum(x, x, x, h, beta, 1:n) - loo * beta[1] * x
    sKy <- ksum(x, y, x, h, beta, 1:n) - loo * beta[1] * y
    sKx2 <- ksum(x, x^2, x, h, beta, 1:n) - loo * beta[1] * x^2
    sKxy <- ksum(x, x * y, x, h, beta, 1:n) - loo * beta[1] * x * y

    # return local-linear estimate of regression function at sample points
    (((sKx2 * sKy - sKx * sKxy) + (sK * sKxy - sKx * sKy) * x) / (sK * sKx2 - sKx^2))
  }
}




# fk_NW computes the local constant (Nadaraya-Watson) regression estimate from sample (x, y)
# Arguments, output, and interpretation of the function is the same as in the case of fk_loc_lin.
# The differences lie only in the fact that the actual regression function estimates have a much
# simpler form (weighted sum of y)/(sum of weights) = sKy/sK


fk_NW <- function(x, y, h, beta, nbin = NULL, ngrid = NULL, from = NULL, to = NULL, loo = 0){
  n <- length(x)
  if(!is.null(nbin)){
    if(is.null(from)) from <- x[1]
    if(is.null(to)) to <- x[n]
    xs <- seq(from, to, length = nbin)
    wts1 <- sm_bin_wts(x, rep(1, n), nbin, xs[1], xs[nbin])
    wtsy <- sm_bin_wts(x, y, nbin, xs[1], xs[nbin])
    sK <- ksum(xs, wts1, xs, h, beta, 1:nbin) - loo * beta[1]
    sK[sK < 1e-20] <- 1e-20
    sKy <- ksum(xs, wtsy, xs, h, beta, 1:nbin) - loo * beta[1] * wtsy / wts1
    sKy / sK
  }
  else if(!is.null(ngrid)){
    if(is.null(from)) from <- x[1]
    if(is.null(to)) to <- x[n]
    xs <- seq(from, to, length = ngrid)
    sK <- ksum(x, rep(1, n), xs, h, beta)
    sKy <- ksum(x, y, xs, h, beta)
    sKy / sK
  }
  else{
    sK <- ksum(x, rep(1, n), x, h, beta, 1:n) - loo * beta[1]
    sK[sK < 1e-20] <- 1e-20
    sKy <- ksum(x, y, x, h, beta, 1:n) - loo * beta[1] * y
    sKy / sK
  }
}


### Methods for class fk_regression

plot.fk_regression <- function(x, main = NULL, ...){
  if(is.null(main)){
    main <- paste('fk_regression(', names(x$call)[2], ' = ', as.character(x$call)[2], sep = '')
    if(length(x$call)>2) for(argnum in 3:length(x$call)) main <- paste(main, ', ', names(x$call)[argnum], ' = ', as.character(x$call)[argnum], sep = '')
    main <- paste(main, ')', sep = '')
  }
  if(length(x$x) > 5000){
    message('Large sample. Only showing a random subset of size 5000')
    sm <- sample(1:length(x$x), 5000)
    plot(x$x[sm], x$y[sm], main = main, xlab = 'x', ylab = 'y', ...)
    if(length(x$x_eval) > 5000){
      sme <- order(sample(1:length(x$x_eval)))
      lines(x$x_eval[sme], x$y_fitted[sme], col = 2, lwd = 2)
    }
    else lines(x$x_eval, x$y_fitted, col = 2, lwd = 2)
  }
  else{
    plot(x$x, x$y, main = main, xlab = 'x', ylab = 'y', ...)
    lines(x$x_eval, x$y_fitted, col = 2, lwd = 2)
  }
}

print.fk_regression <- function(x, ...){
  cat('Call: \n \n')
  print(x$call)
  cat('\n')
  cat(paste('Data: (x, y) (', length(x$x), ' obs.); bandwidth = ', round(x$h, 4), ' \n \n', sep = ''))
}


predict.fk_regression <- function(object, xtest = NULL, ...){
  if(is.null(xtest)) xtest <- object$x

  mn <- mean(object$x)

  xtest <- xtest - mn
  x <- object$x - mn

  o <- order(x)

  otest <- order(xtest)

  if(object$type == 'loc-lin'){
    sK <- ksum(x[o], numeric(length(x)) + 1, xtest[otest], object$h, object$betas)
    sKx <- ksum(x[o], x[o], xtest[otest], object$h, object$betas)
    sKy <- ksum(x[o], object$y[o], xtest[otest], object$h, object$betas)
    sKx2 <- ksum(x[o], x[o]^2, xtest[otest], object$h, object$betas)
    sKxy <- ksum(x[o], x[o] * object$y[o], xtest[otest], object$h, object$betas)
    yhat <- (((sKx2 * sKy - sKx * sKxy) + (sK * sKxy - sKx * sKy) * xtest[otest]) / (sK * sKx2 - sKx^2))[rank(xtest)]
  }
  else{
    sK <- ksum(x[o], numeric(length(x)) + 1, xtest[otest], object$h, object$betas)
    sKy <- ksum(x[o], object$y[o], xtest[otest], object$h, object$betas)
    yhat <- (sKy / sK)[rank(xtest)]
  }
  yhat
}
