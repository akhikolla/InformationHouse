#' Geometrically Designed Spline regression estimation
#'
#' @description \code{NGeDS} constructs a Geometrically Designed  variable
#' knots spline regression model  referred to as
#' a GeDS model, for a response having a Normal distribution.
#'
#' @export
#'
#' @param formula a description of the structure of the model to be fitted, including the
#' dependent and independent variables. See \code{\link[=formula.GeDS]{formula}} for details.
#' @param data an optional data frame, list or environment containing the variables of the model.
#' If not found in \code{data}, the variables are taken from \code{environment(formula)},
#' typically the environment from which \code{NGeDS} is called.
#' @param weights an optional vector of `prior weights' to be put on
#' the observations in the fitting
#' process in case the user requires weighted GeDS fitting.
#' It should be \code{NULL} or a numeric vector of the same length as the response
#'  variable in the argument \code{\link[=formula.GeDS]{formula}}.
#' @param beta numeric parameter in the interval \eqn{[0,1]}
#' tuning the knot placement in stage A of GeDS. See details.
#' @param phi numeric parameter in the interval \eqn{[0,1]} specifying the threshold for
#'  the stopping rule  (model selector) in stage A of GeDS. See also \code{stoptype} and details below.
#' @param min.intknots optional parameter allowing the user to set
#' a minimum number of internal knots required. By default equal to zero.
#' @param max.intknots optional parameter allowing the user to set a maximum number
#'  of internal knots to be added by the GeDS estimation algorithm.
#' By default equal to the number of knots for the saturated GeDS model.
#' @param q numeric parameter which allows to fine-tune the stopping rule of stage A of GeDS, by default equal to 2.
#' See details.
#' @param Xextr numeric vector of 2 elements representing the left-most and right-most limits
#' of the interval embedding the observations of the first independent variable. See details.
#' @param Yextr numeric vector of 2 elements representing the left-most and right-most limits
#' of the interval embedding the observations of the second independent variable
#' (if the bivariate GeDS is run). See details.
#' @param show.iters logical variable indicating whether or not to print
#' information at each step.
#' @param stoptype a character string indicating the type of GeDS stopping rule to
#' be used. It should be either one of \code{"SR"}, \code{"RD"} or
#' \code{"LR"}, partial match allowed. See details.
#' @return \code{\link{GeDS-Class}} object, i.e. a list of items that summarizes the main
#' details of the fitted GeDS regression. See \code{\link{GeDS-Class}} for details.
#' Some S3 methods are available in order to make these objects tractable, such as \code{\link[=coef.GeDS]{coef}},
#' \code{\link[=deviance.GeDS]{deviance}}, \code{\link[=knots.GeDS]{knots}}, \code{\link[=predict.GeDS]{predict}} and \code{\link[=print.GeDS]{print}}
#' as well as S4 methods for \code{\link[=lines.GeDS]{lines}} and \code{\link[=plot.GeDS]{plot}}.
#' @details
#' The  \code{NGeDS} function implements the GeDS methodology, recently developed by Kaishev et al. (2016)
#'  and extended
#' in the \code{\link{GGeDS}} function for the more general GNM, (GLM) context, allowing for the response to have any
#' distribution from the Exponential Family.
#' Under the GeDS approach the (non-)linear predictor is viewed as a spline with
#' variable knots which are estimated along with the regression coefficients and the order of the spline, using a
#' two stage algorithm.
#' In stage A, a linear variable-knot spline is fitted to the data applying iteratively
#' least squares  regression (see \code{\link[stats]{lm}} function). In stage B, a Schoenberg variation diminishing spline approximation to the fit
#' from stage A is constructed, thus simultaneously producing spline fits of order 2,
#' 3 and 4, all of which are included in the output, a \code{\link{GeDS-Class}} object.
#'
#' As noted in \code{\link[=formula.GeDS]{formula}}, the argument \code{formula} allows the user to specify models
#' with two components, a
#'  spline regression (non-parametric) component involving part of the independent variables identified through
#'   the function \code{f}
#'  and an optional  parametric component involving the remaining independent variables.
#'  For \code{NGeDS} one or two independent variables are allowed for the spline component
#'  and arbitrary many independent variables for the parametric component.
#' Failure to specify the independent variable for the  spline regression component through
#' the function \code{f} will return an error.
#'  See \code{\link[=formula.GeDS]{formula}}.
#'
#'  Within the argument \code{formula}, similarly as in other R functions, it is possible to
#'  specify one or more offset variables, i.e. known terms with fixed regression coefficients equal to 1.
#'  These terms should be identified via the function \code{\link[stats]{offset}}.
#'
#'  The parameter \code{beta} tunes the placement of a new knot in stage A of the algorithm.
#'  Once a current second-order  spline is fitted to the data the
#'  regression residuals are computed and grouped by their sign.
#'  A new knot is placed  at a location defined by the group
#'  for which a certain measure attains its maximum.
#'   The latter measure is   defined  as a weighted linear combination of the range
#'  of each group and the  mean of the absolute residuals within it.
#'  The parameter \code{beta} determines the weights in this measure correspondingly as  \code{beta}
#'  and \code{1 - beta}.
#'  The  higher it is, the more weight is put to the mean of the residuals and the less to the range
#'  of their corresponding x-values. The default value of \code{beta} is 0.5.
#'
#'  The argument \code{stoptype} allows to choose between three alternative stopping rules for the knot selection
#'  in stage A of GeDS, the \code{"RD"}, that stands for \emph{Ratio of Deviances}, the \code{"SR"},
#'  that stands for \emph{Smoothed Ratio} of deviances and the \code{"LR"}, that stands for \emph{Likelihood Ratio}.
#'  The latter is based on the difference of deviances rather than on their ratio as in the case of
#'  \code{"RD"} and \code{"SR"}. Therefore \code{"LR"} can be viewed as
#'  a log likelihood ratio test performed at each iteration of the knot placement.
#'  In each of these cases the corresponding stopping criterion is compared with a threshold value
#'  \code{phi} (see below).
#'
#' The argument \code{phi} provides a threshold value required for the stopping rule to exit
#' the knot placement in stage A of GeDS.
#' The higher the value of \code{phi}, the more knots are added under the \code{"RD"} and \code{"SR"} stopping rules
#' contrary to the case of the stopping rule \code{"LR"} where the lower \code{phi} is,
#' more knots are included in the spline regression. Further details for
#' each of the three alternative stopping rules can be found in Dimitrova et al. (2017).
#'
#'  %In all the three cases, the input parameter \code{q} controls the computation of the ratio
#'  %or difference of deviances. The new deviance is compared with the one computed \code{q} iterations
#'  %before,
#'
#'
#'   The argument \code{q} is an input parameter that allows to fine-tune the stopping rule in stage A.
#'    It identifies the number of consecutive iterations
#'   over which the deviance should exhibit stable convergence so as the knot placement in stage A is terminated.
#'   More precisely,
#'   under any of the rules \code{"RD"}, \code{"SR"} or \code{"LR"}
#'   the deviance at the current iteration is compared to the deviance computed
#'    \code{q} iterations before,  i.e. before selecting
#'   the last \code{q} knots. Setting a higher \code{q} will lead to more knots
#'   being added before exiting stage A of GeDS.
#'
#'
#'
#' @references
#' Kaishev, V.K., Dimitrova, D.S., Haberman, S. and Verrall, R.J. (2016).
#' Geometrically designed, variable knot regression splines.
#' \emph{Computational Statistics}, \strong{31}, 1079--1105. \cr
#' DOI: \href{https://doi.org/10.1007/s00180-015-0621-7}{doi.org/10.1007/s00180-015-0621-7}
#'
#' Dimitrova, D.S., Kaishev, V.K., Lattuada A. and Verrall, R.J. (2017).
#' Geometrically designed, variable knot splines in Generalized (Non-)Linear Models
#'
#'
#' @seealso \link{GGeDS}; \link{GeDS-Class}; S3 methods such as \link{coef.GeDS}, \link{deviance.GeDS}, \link{knots.GeDS},
#' \link{print.GeDS} and \link{predict.GeDS}; \link{Integrate} and \link{Derive}; \link{PPolyRep}.
#'
#' @examples
#' ###################################################
#' # Generate a data sample for the response variable
#' # Y and the single covariate X
#' set.seed(123)
#' N <- 500
#' f_1 <- function(x) (10*x/(1+100*x^2))*4+4
#' X <- sort(runif(N, min = -2, max = 2))
#' # Specify a model for the mean of Y to include only a component
#' # non-linear in X, defined by the function f_1
#' means <- f_1(X)
#' # Add (Normal) noise to the mean of Y
#' Y <- rnorm(N, means, sd = 0.1)
#'
#' # Fit a Normal GeDS regression using NGeDS
#' (Gmod <- NGeDS(Y ~ f(X), beta = 0.6, phi = 0.995, Xextr = c(-2,2)))
#'
#' # Apply some of the available methods, e.g.
#' # coefficients, knots and deviance extractions for the
#' # quadratic GeDS fit
#' # Note that the first call to the function knots returns
#' # also the left and right limits of the interval containing
#' # the data
#' coef(Gmod, n = 3)
#' knots(Gmod, n = 3)
#' knots(Gmod, n = 3, options = "internal")
#' deviance(Gmod, n = 3)
#'
#' # Add a covariate, Z, that enters linearly
#' Z <- runif(N)
#' Y2 <- Y + 2*Z + 1
#' # Re-fit the data using NGeDS
#' (Gmod2 <- NGeDS(Y2 ~ f(X) + Z, beta = 0.6, phi = 0.995, Xextr = c(-2,2)))
#' coef(Gmod2, n = 3)
#' coef(Gmod2, onlySpline = FALSE, n = 3)
#'
#' \dontrun{
#' ##########################################
#' # Real data example
#' # See Kaishev et al. (2016), section 4.2
#' data('BaFe2As2')
#' (Gmod2 <- NGeDS(intensity ~ f(angle), data = BaFe2As2, beta = 0.6, phi = 0.99, q = 3))
#' plot(Gmod2)
#' }
#'
#' #########################################
#' # bivariate example
#' # See Dimitrova et al. (2017), section 5
#'
#' # Generate a data sample for the response variable
#' # Z and the covariates X and Y assuming Normal noise
#' set.seed(123)
#' doublesin <- function(x){
#'  sin(2*x[,1])*sin(2*x[,2])
#' }
#'
#' x <- (round(runif(400, min = 0, max = 3),2))
#' y <- (round(runif(400, min = 0, max = 3),2))
#' z <- doublesin(cbind(x,y))
#' z <- z+rnorm(400, 0, sd = 0.1)
#' # Fit a two dimensional GeDS model using NGeDS
#' (BivGeDS <- NGeDS(z ~ f(x, y) , phi = 0.9, beta = 0.3,
#' Xextr = c(0, 3), Yextr = c(0, 3)))
#'
NGeDS<- function(formula, data,   weights,
                beta=.5, phi = 0.99, min.intknots = 0,
                max.intknots = 500, q = 2, Xextr=NULL,
                Yextr=NULL, show.iters=FALSE,stoptype="RD"){
  save <- match.call()
  if (missing(data))    data <- environment(formula)
  wn <- match("weights",names(save),0L)
  wn <- save[c(1L,wn)]
  wn[[1L]] <- quote(list)
  weights <- eval(wn,data,parent.frame())$weights

  newdata <- read.formula(formula, data)
  X <- newdata$X
  Y <- newdata$Y
  offset <- newdata$offset
  Z <- newdata$Z
  if(!is.null(Z)) Z <- as.matrix(Z)

  if(is.null(weights)) weights=rep(1,NROW(X))
  weights <- as.numeric(weights)
  if(is.null(Z)) {
    if(NROW(X)!=length(Y) || NROW(X)!=length(weights)) stop("length of 'X', 'Y' and 'weights' must match")
  } else {
    if(NROW(X)!=length(Y) || NROW(X)!=NROW(Z) || NROW(X)!=length(weights)) stop("length of 'X', 'Y', 'Z' and 'weights' must match")
  }
  if(is.null(X) || is.null(Y)) stop("Null arguments cannot be passed to this function")
  if(length(phi) != 1 || length(beta) != 1 || length(max.intknots) != 1 ||
     length(q) != 1 || length(show.iters) != 1)
    stop("'phi', 'beta',  'max.intknots', 'q' and 'show.iters' must have length 1")
  if(beta > 1 || beta < 0) stop("'beta' should be a value in [0,1]")
    phi <- as.numeric(phi)
    if(phi > 1 || phi < 0) stop("'phi' should be a value in [0,1]")
  RSSnew <- numeric()
  max.intknots <- as.integer(max.intknots)
  min.intknots <- as.integer(min.intknots)
  q <- as.integer(q)
  show.iters <- as.logical(show.iters)
  if(!is.null(Xextr)){ if(length(Xextr) !=2) stop("'Xextr' must have length 2")}
  if(!is.null(Yextr)){ if(length(Yextr) !=2) stop("'Yextr' must have length 2")}

  if(ncol(X)==1){
    if(is.null(Z)){
      if(anyNA(X) || anyNA(Y) || anyNA(weights) || anyNA(offset)) {
        warning("NAs deleted from 'X', 'Y', 'offset' and 'weights'")
        tmp <- is.na(X) | is.na(Y) | is.na(weights) | is.na(offset)
        Y <- Y[!tmp]
        X <- X[!tmp]
        weights <- weights[!tmp]
        offset <- offset[!tmp]
      }
      idx <- order(X)
      X <- X[idx]
      Y <- Y[idx]
      weights <- weights[idx]
      offset <- offset[idx]
    } else {
      if(anyNA(X) || anyNA(Y) || anyNA(weights) || anyNA(Z) || anyNA(offset)) {
        warning("NAs deleted from 'X', 'Y', 'Z', 'offset' and 'weights'")
        tmp <- apply(Z, anyNA,MARGIN=1)
        tmp <- is.na(X) | is.na(Y) | is.na(weights) | is.na(offset) | tmp
        Y <- Y[!tmp]
        X <- X[!tmp]
        weights <- weights[!tmp]
        Z <- Z[!tmp,]
        offset <- offset[!tmp]
      }
      idx <- order(X)
      X <- X[idx]
      Y <- Y[idx]
      Z <- Z[idx,1:NCOL(Z)]
      weights <- weights[idx]
      offset <- offset[idx]
    }
    if (any(weights < 0))  stop("Negative weights not allowed")
    Xextr <-if (is.null(Xextr)) range(X) else as.numeric(Xextr)
    out <- UnivariateFitter(X=X, Y=Y, Z=Z, offset=offset,
                                     weights = weights,
                                     beta=beta, phi = phi,
                                     min.intknots = min.intknots,
                                     max.intknots = max.intknots, q = q,
                                     extr=Xextr, show.iters=show.iters,stoptype=stoptype)
  } else {
    if(ncol(X)==2){
      if(is.null(Z)){
        if(anyNA(X) || anyNA(Y) || anyNA(weights)) {
          warning("NAs deleted from 'X', Y' and 'weights'")
          tmp <- apply(X, anyNA,MARGIN=1)
          tmp <- tmp | is.na(Y) | is.na(weights)
          Y[tmp] <- NULL
          X[tmp,] <- X[!tmp,]
          weights <- weights[!tmp]
        }
      } else {
        if(anyNA(X) || anyNA(Y) || anyNA(weights) || anyNA(Z)) {
          warning("NAs deleted from 'X', 'Y', 'Z' and 'weights'")
          tmp <- apply(Z, anyNA,MARGIN=1)
          tmp1 <- apply(X, anyNA,MARGIN=1)
          tmp <- tmp1 | is.na(Y) | is.na(weights) | tmp
          Y <- Y[!tmp]
          X <- X[!tmp,]
          Z <- Z[!tmp,]
          weights <- weights[!tmp]
        }
      }
      if (any(weights < 0))  stop("Negative weights not allowed")
      Indicator <- table(X[,1],X[,2])
      Yextr <-if (is.null(Yextr)) range(X[,2]) else as.numeric(Yextr)
      Xextr <-if (is.null(Xextr)) range(X[,1]) else as.numeric(Xextr)
      out <- BivariateFitter( Y=X[,2], X=X[,1],W=Z, Z=Y,  weights = weights,
                                       Indicator=Indicator, beta=beta, alpha = phi,
                                       min.intknots = min.intknots,
                                       max.intknots = max.intknots, q = q,
                                       Xextr=Xextr,
                                       Yextr=Yextr, show.iters=show.iters)
    } else {
      stop("Incorrect number of columns of the independent variable")
    }
  }
  out$Formula <- formula
  out$extcall <- save
  out$terms <- newdata$terms
  out$znames <- getZnames(newdata)
  return(out)
}

