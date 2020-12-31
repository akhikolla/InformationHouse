#' @export
#' @title IRLS Estimation
#' @param x a matrix of regression functions (e.g. B-splines and/or terms of the parametric part)
#' evaluated at the sample values of the covariate(s).
#' @param y a vector of size \eqn{N} containing the observed values of the response variable \eqn{y}.
#' @param weights an optional vector of `prior weights' to be put on
#' the observations in case the user requires weighted IRLS fitting.
#'  It is a vector of 1s by default.
#' @param mustart initial values for the vector of means of the response variable in the IRLS regression
#' estimation. Must be a vector of length \eqn{N}.
#' @param offset a vector of size \eqn{N} that can be used to specify a fixed covariate
#' to be included in the predictor model  avoiding the estimation of its corresponding regression coefficient.
#' In case  more than one covariate is fixed, the user should sum the corresponding coordinates of the fixed covariates
#'  to produce one common \eqn{N}-vector of coordinates.
#' @param family a description of the error distribution and link function to be used in the model. This can be a
#' character string naming a family function (e.g. \code{"gaussian"}),
#' the family function itself (e.g. \code{\link[stats]{gaussian}})
#' or the result of a call to a family function (e.g. \code{gaussian()}).
#'  See \link[stats]{family} for details on family functions.
#' @param control a list of parameters for controlling the IRLS fitting process to be passed
#' on to \code{\link[stats]{glm.control}}. See \code{\link[stats]{glm.fit}} for further details.
#' @details
#' This function is a slightly modified version of the \code{\link[stats]{glm.fit}} from the package \pkg{stats} to which
#' we refer for further details.
#' The difference in the inputs of \code{IRLSfit} and \code{\link[stats]{glm.fit}} is that the former admits initial values
#' only for the vector of means.
#'
#' The output from \code{IRLSfit} has some additional slots compared to \code{\link[stats]{glm.fit}}.
#' We note that the slots \code{weights}, \code{res2} and \code{z} contain
#' values of the IRLS weights, ``working residuals" and transformed responses computed \emph{after} the last IRLS iteration, i.e.
#'  they are based on the estimated coefficients that are returned by \code{IRLSfit}.
#'
#' The source code of \code{IRLSfit} contains also some commented lines that produce useful plots at each IRLS iteration.
#' Normally, printing these plots is time consuming, but they could be run for inspection purposes.
#'
#'
#' @description This function is an implementation of the IRLS estimation algorithm adjusted to the specific usage in the  function
#' \code{\link{SplineReg_GLM}}.
#'
#'
#' @return A list containing:
#' \item{coefficients}{a named vector containing the estimated regression coefficients;}
#' \item{residuals}{the `working' residuals, that are the residuals in the final iteration of the IRLS fit. Since cases with
#' zero weights are omitted, their working residuals are \code{NA};}
#' \item{res2}{the working residuals after the final IRLS iteration. They are used within the  knot placement steps
#' of stage A of GeDS;}
#' \item{fitted.values}{the fitted mean values, obtained by transforming the predictor by the inverse of the link function;}
#' \item{rank}{the numeric rank of the fitted linear model;}
#' \item{family}{the \code{\link[stats]{family}} object used;}
#' \item{linear.predictors}{the fitted predictor;}
#' \item{deviance}{a vector containing the deviances obtained at each IRLS iteration;}
#' \item{lastdeviance}{the deviance at the last IRLS iteration;}
#' \item{null.deviance}{The deviance for the null model (see \code{\link[stats]{glm}} documentation);}
#' \item{iter}{the number of IRLS iterations performed;}
#' \item{weights}{the working weights after the last IRLS iteration;}
#' \item{prior.weights}{the ``prior weights" (see the \code{weights} argument);}
#' \item{df.residual}{the residual degrees of freedom;}
#' \item{df.null}{the residual degrees of freedom for the null model;}
#' \item{y}{the vector of values of the response variable used in the fitting;}
#' \item{z}{the transformed responses computed after the last IRLS iteration;}
#' \item{converged}{logical. Was the IRLS algorithm judged to have converged?}
#' \item{boundary}{logical. Is the fitted value on the boundary of the attainable values?}
#' In addition, non-empty fits will have components \code{qr}, \code{R} and \code{effects} relating to the final weighted
#' linear fit, see \code{\link{lm.fit}} documentation.
#'
#' @seealso \code{\link[stats]{glm.fit}}
IRLSfit <- function (x, y, weights = rep(1, nobs),
                     mustart = NULL, offset = rep(0, nobs), family = gaussian(),
                     control = list())
{
  start = NULL
  devi2 <- NULL
  flag <- FALSE
  control <- do.call("glm.control", control)
  x <- as.matrix(x)
  xnames <- dimnames(x)[[2L]]
  ynames <- if (is.matrix(y))
    rownames(y)
  else names(y)
  conv <- FALSE
  nobs <- NROW(y)
  #  if(is.null(offset)) offset=rep(0, nobs)
  n <- NULL

  nvars <- ncol(x)
  EMPTY <- nvars == 0
  if (is.null(weights))
    weights <- rep.int(1, nobs)
  if (is.null(offset))
    offset <- rep.int(0, nobs)
  variance <- family$variance
  linkinv <- family$linkinv
  if (!is.function(variance) || !is.function(linkinv))
    stop("'family' argument seems not to be a valid family object",
         call. = FALSE)
  dev.resids <- family$dev.resids
  aic <- family$aic
  mu.eta <- family$mu.eta



  unless.null <- function(x, if.null) if (is.null(x))
    if.null else x
  valideta <- unless.null(family$valideta, function(eta) TRUE)
  validmu <- unless.null(family$validmu, function(mu) TRUE)
  if (is.null(mustart)) {
    eval(family$initialize)
  }  else {
    mukeep <- mustart

    eval(family$initialize)

    mustart <- mukeep
  }
  if (EMPTY) {
    eta <- rep.int(0, nobs) + offset
    if (!valideta(eta))
      stop("invalid linear predictor values in empty model",
           call. = FALSE)
    mu <- linkinv(eta)
    if (!validmu(mu))
      stop("invalid fitted means in empty model", call. = FALSE)
    dev <- sum(dev.resids(y, mu, weights))
    w <- ((weights * mu.eta(eta)^2)/variance(mu))^0.5
    residuals <- (y - mu)/mu.eta(eta)
    good <- rep_len(TRUE, length(residuals))
    boundary <- conv <- TRUE
    coef <- numeric()
    iter <- 0L
  }
  else {
    coefold <- NULL
    eta <- family$linkfun(mustart)
    mu <- mustart
    if (!(validmu(mu) && valideta(eta)))
      stop("cannot find valid starting values: please specify some",
           call. = FALSE)
    devold <- sum(dev.resids(y, mu, weights))
    boundary <- conv <- FALSE

    devi2 <- devold

    for (iter in 1L:control$maxit) {
      good <- weights > 0
      varmu <- variance(mu)[good]
      if (anyNA(varmu))
        stop("NAs in V(mu)")
      if (any(varmu == 0))
        stop("0s in V(mu)")
      mu.eta.val <- mu.eta(eta) #############################
      if (any(is.na(mu.eta.val[good])))
        stop("NAs in d(mu)/d(eta)")








      good <- (weights > 0) & (mu.eta.val != 0)

      if (all(!good)) {
        conv <- FALSE
        warning(gettextf("no observations informative at iteration %d",
                         iter), domain = NA)
        break
      }

      z <- (eta - offset)[good] + (y - mu)[good]/mu.eta.val[good]


      #zeroes <- !(z>0)
      #z[zeroes] <- 0
      w <- sqrt((weights[good] * mu.eta.val[good]/variance(mu)[good])*mu.eta.val[good])

      ## plots whithin IRLS uncommenting following lines
      #if(iter %in% c(1L ,2L ,3L)) {
      #  if(iter==1L) {
      #    par(mfrow=c(2,2),mai=c(0.5,0.5,0.5,0.5))}
      #  buoni <- !is.infinite(family$linkfun(y))
      #  ys <- family$linkfun(y[buoni])
      #  rangey <- if(any(!buoni)) c(-max(abs(ys)),max(abs(ys))) else NULL
      #  plot(X[buoni],ys,xlim=range(X),ylim=rangey,col="red")
      #  points(X[good],z)
      #  lines(X[good],eta,col="blue")
      #
      #}
      ## other plots
      #if(iter==3L)  {       plot(X[buoni],ys,xlim=range(X),ylim=rangey,col="red")
      #points(X[good],(eta - offset)[good] + first.term)
      #lines(X[good],eta,col="blue")
      #dev.off()
      #}
      #      print(iter)
      fit <- .lm.fit(x[good, , drop = FALSE] *
                       w, z * w, min(1e-07, control$epsilon/1000))


      if (any(!is.finite(fit$coefficients))) {
        conv <- FALSE
        warning(gettextf("non-finite coefficients at iteration %d",
                         iter), domain = NA)
        break
      }
      if (nobs < fit$rank)
        stop(sprintf(ngettext(nobs, "X matrix has rank %d, but only %d observation",
                              "X matrix has rank %d, but only %d observations"),
                     fit$rank, nobs), domain = NA)
      start[fit$pivot] <- fit$coefficients
      eta <- drop(x %*% start)

      mu <- linkinv(eta <- eta + offset)

      ## nice plots uncommenting following lines
      #if(T){
      ##  Sys.sleep(.05)  ## will allow to inspect the plots, but the code becomes awfully slow


      #par(mfrow=c(2,1),mai=c(0.5,0.5,0.5,0.5))
      #print(round(start,3))

      #plot(X,y)
      #lines(X[good],mu,col="red")
      #buoni <- !is.infinite(family$linkfun(y))
      #ys <- family$linkfun(y[buoni])
      #rangey <- if(any(!buoni)) c(-max(abs(ys)),max(abs(ys))) else NULL
      #plot(X[buoni],ys,xlim=range(X),ylim=rangey,col="red")
      #points(X[good],z)
      #lines(X[good],eta,col="red")
      #}



      dev <- sum(dev.resids(y, mu, weights))
      devi2 <-  c(devi2,dev)
      if (control$trace)
        cat("Deviance = ", dev, " Iterations - ", iter,
            "\n", sep = "")
      boundary <- FALSE
      if (!is.finite(dev)) {
        if (is.null(coefold)){
          if(flag){
            stop("no valid set of coefficients has been found: please supply starting values",
                 call. = FALSE)
          } else {
            flag <- TRUE
            eval(family$initialize)

            eta <- family$linkfun(mustart)
            mu <- mustart
            next
          }
          } else {flag <- FALSE}
        warning("step size truncated due to divergence",
                call. = FALSE)
        ii <- 1
        while (!is.finite(dev)) {
          if (ii > 100 )#control$maxit)
            stop("inner loop 1; cannot correct step size",
                 call. = FALSE)
          ii <- ii + 1
          start <- (start + coefold)/2
          eta <- drop(x %*% start)
          mu <- linkinv(eta <- eta + offset)
          dev <- sum(dev.resids(y, mu, weights))
        }
        boundary <- TRUE
        if (control$trace)
          cat("Step halved: new deviance = ", dev, "\n",
              sep = "")
      }
      if (!(valideta(eta) && validmu(mu))) {
        if (is.null(coefold)){
          if(flag){
            stop("no valid set of coefficients has been found: please supply starting values",
                 call. = FALSE)
          } else {
            flag <- TRUE
            eval(family$initialize)
            eta <- family$linkfun(mustart)
            mu <- mustart
            next
          }
        } else {flag <- FALSE}
        warning("step size truncated: out of bounds",
                call. = FALSE)
        ii <- 1
        while (!(valideta(eta) && validmu(mu))) {
          if (ii > control$maxit)
            stop("inner loop 2; cannot correct step size",
                 call. = FALSE)
          ii <- ii + 1
          start <- (start + coefold)/2
          eta <- drop(x %*% start)
          mu <- linkinv(eta <- eta + offset)
        }
        boundary <- TRUE
        dev <- sum(dev.resids(y, mu, weights))
        if (control$trace)
          cat("Step halved: new deviance = ", dev, "\n",
              sep = "")
      }
      if (abs(dev - devold)/(0.1 + abs(dev)) < control$epsilon) {
        conv <- TRUE
        coef <- start
        break
      }
      else {
        devold <- dev
        coef <- coefold <- start
      }
    }
    if (!conv)
      warning("IRLS algorithm did not converge", call. = FALSE)
    if (boundary)
      warning("IRLS algorithm stopped at boundary value",
              call. = FALSE)
    eps <- 10 * .Machine$double.eps
    if (family$family == "binomial") {
      if (any(mu > 1 - eps) || any(mu < eps))
        warning("IRLS fitted probabilities numerically 0 or 1 occurred",
                call. = FALSE)
    }
    if (family$family == "poisson") {
      if (any(mu < eps))
        warning("IRLS fitted rates numerically 0 occurred",
                call. = FALSE)
    }
    if (fit$rank < nvars)
      coef[fit$pivot][seq.int(fit$rank + 1, nvars)] <- NA
    xxnames <- xnames[fit$pivot]
    residuals <- (y - mu)/mu.eta(eta)
    fit$qr <- as.matrix(fit$qr)
    nr <- min(sum(good), nvars)
    if (nr < nvars) {
      Rmat <- diag(nvars)
      Rmat[1L:nr, 1L:nvars] <- fit$qr[1L:nr, 1L:nvars]
    }
    else Rmat <- fit$qr[1L:nvars, 1L:nvars]
    Rmat <- as.matrix(Rmat)
    Rmat[row(Rmat) > col(Rmat)] <- 0
    names(coef) <- xnames
    colnames(fit$qr) <- xxnames
    dimnames(Rmat) <- list(xxnames, xxnames)
  }
  names(residuals) <- ynames
  names(mu) <- ynames
  names(eta) <- ynames
  wt <- rep.int(0, nobs)
  wt[good] <- w^2
  names(wt) <- ynames
  names(weights) <- ynames
  names(y) <- ynames
  if (!EMPTY)
    names(fit$effects) <- c(xxnames[seq_len(fit$rank)], rep.int("",
                                                                sum(good) - fit$rank))
  wtdmu <- sum(weights * y)/sum(weights)
  nulldev <- sum(dev.resids(y, wtdmu, weights))
  n.ok <- nobs - sum(weights == 0)
  nulldf <- n.ok - 1
  rank <- if (EMPTY)
    0
  else fit$rank
  resdf <- n.ok - rank


  good <- weights > 0
  varmu <- variance(mu)[good]
  if (anyNA(varmu))
    stop("NAs in V(mu)")
  if (any(varmu == 0))
    stop("0s in V(mu)")
  mu.eta.val <- mu.eta(eta)
  if (any(is.na(mu.eta.val[good])))
    stop("NAs in d(mu)/d(eta)")
  good <- (weights > 0) & (mu.eta.val != 0)

  z <-(eta-offset)[good] + (y - mu)[good]/mu.eta.val[good] #

  res2 <- (y - mu)[good]/mu.eta.val[good]

  w <- ((mu.eta(eta)^2)/variance(mu))^0.5
  wt <- rep.int(0, nobs)
  wt[good] <- w^2

  if(!conv) print("IRLS did not converge.")
  list(coefficients = coef, residuals = residuals, res2 = res2, fitted.values = mu,
       effects = if (!EMPTY) fit$effects, R = if (!EMPTY) Rmat, rank = rank,
       qr = if (!EMPTY) structure(fit[c("qr", "rank","qraux", "pivot", "tol")],class = "qr"),
       family = family, linear.predictors = eta, deviance = devi2,
       lastdeviance = dev, null.deviance = nulldev, iter = iter, weights = wt, prior.weights = weights,
       df.residual = resdf, df.null = nulldf, y = y, z = z, converged = conv,
       boundary = boundary)
}



