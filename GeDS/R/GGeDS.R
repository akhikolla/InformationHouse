

#' Generalized Geometrically Designed Spline regression estimation
#'
#' @description \code{GGeDS} constructs a Geometrically Designed (univariate) variable
#' knots spline regression model for the predictor  in the context of Generalized (Non-)Linear Models,
#' referred to as a GeDS model, for a response with a pre-specified distribution from
#' the Exponential Family.
#'
#' @export
#'
#' @param formula a description of the structure of the predictor model to be fitted, including the
#' dependent and independent variables. See \code{\link[=formula.GeDS]{formula}} for details.
#' @param data an optional data frame, list or environment containing the variables of the predictor model.
#' In case the variables are not found in \code{data}, they are taken from \code{environment(formula)},
#' typically the environment from which \code{GGeDS} is called.
#' @param family a description of the error distribution and link function to be used in the model. This can be a
#' character string naming a family function (e.g. \code{"gaussian"}),
#' the family function itself (e.g. \code{\link[stats]{gaussian}})
#' or the result of a call to a family function (e.g. \code{gaussian()}).
#'  See \link[stats]{family} for details on family functions.
#' @param weights an optional vector of `prior weights' to be put on
#' the observations in the fitting
#' process in case the user requires weighted GeDS fitting.
#' It is \code{NULL} by default.
#' @param beta numeric parameter in the interval \eqn{[0,1]}
#' tuning the knot placement in stage A of GeDS. See details below.
#' @param phi numeric parameter in the interval \eqn{[0,1]} specifying the threshold for
#'  the stopping rule  (model selector) in stage A of GeDS. See also \code{stoptype} and details below.
#' @param min.intknots optional parameter allowing the user to set
#' a minimum number of internal knots required. By default equal to zero.
#' @param max.intknots optional parameter allowing the user to set a maximum number
#'  of internal knots to be added by the GeDS estimation algorithm.
#' By default equal to the number of knots for the saturated GeDS model (i.e. \eqn{N-2},
#' where \eqn{N} is the number of observations).
#' @param q numeric parameter which allows to fine-tune the stopping rule of stage A of GeDS, by default equal to 2.
#' See details below.
#' @param Xextr numeric vector of 2 elements representing the left-most and right-most limits
#' of the interval embedding the observations of the independent variable. See details.
#' @param show.iters logical variable indicating whether or not to print
#' information at each step. By default equal to \code{FALSE}.
#' @param stoptype a character string indicating the type of GeDS stopping rule to
#' be used. It should be either one of \code{"SR"}, \code{"RD"} or
#' \code{"LR"}, partial match allowed. See details below.
#' @return A \code{\link{GeDS-Class}} object, i.e. a list of items that summarizes  the main
#' details of the fitted GeDS regression. See \code{\link{GeDS-Class}} for details.
#' Some S3 methods are available in order to make these objects tractable, such as \code{\link[=coef.GeDS]{coef}},
#' \code{\link[=deviance.GeDS]{deviance}}, \code{\link[=knots.GeDS]{knots}}, \code{\link[=predict.GeDS]{predict}} and \code{\link[=print.GeDS]{print}}
#' as well as S4 methods for \code{\link[=lines.GeDS]{lines}} and \code{\link[=plot.GeDS]{plot}}.
#' @details
#' The  \code{GGeDS} function extends the GeDS methodology, recently developed by Kaishev et al. (2016) and implemented
#' in the \code{\link{NGeDS}} function for the Normal case,
#' to the more general GNM (GLM) context, allowing for the response to have any
#' distribution from the Exponential Family.
#' Under the GeDS-GNM approach the (non-)linear predictor is viewed as a spline with
#' variable knots which are estimated along with the regression coefficients and the
#' order of the spline, using a two stage procedure.
#' In stage A, a linear variable-knot spline is fitted to the data applying iteratively
#' re-weighted least squares (see \code{\link{IRLSfit}} function). In stage B, a Schoenberg variation
#'  diminishing spline approximation to the fit
#' from stage A is constructed, thus simultaneously producing spline fits of order 2,
#' 3 and 4, all of which are included in the output, a \code{\link{GeDS-Class}} object.
#' A detailed description of the underlying algorithm can be found in Dimitrova et al. (2017).
#'
#' As noted in \code{\link[=formula.GeDS]{formula}}, the argument \code{formula} allows the user to specify predictor models
#' with two components, a
#'  spline regression (non-parametric) component involving part of the independent variables identified through
#'   the function \code{f}
#'  and an optional  parametric component involving the remaining independent variables.
#'  For \code{GGeDS} only one independent variable is allowed for the spline component and arbitrary many independent variables for the
#'  parametric component of the predictor.
#' Failure to specify the independent variable for the  spline regression component through
#' the function \code{f} will return an error.
#'  See \code{\link[=formula.GeDS]{formula}}.
#'
#'  Within the argument \code{formula}, similarly as in other \R functions, it is possible to
#'  specify one or more offset variables, i.e. known terms with fixed regression coefficients equal to 1.
#'  These terms should be identified via the function \code{\link[stats]{offset}}.
#'
#'  The parameter \code{beta} tunes the placement of a new knot in stage A of the algorithm.
#'  Once a current second-order  spline is fitted to the data the
#'  'working' residuals (see \code{\link{IRLSfit}}) are computed and grouped by their sign.
#'  A new knot is placed  at a location defined by the group
#'  for which a certain measure attains its maximum.
#'   The latter measure is   defined  as a weighted linear combination of the range
#'  of each group and the  mean of the absolute residuals within it.
#'  The parameter \code{beta} determines the weights in this measure correspondingly as  \code{beta}
#'  and \code{1 - beta}.
#'  The  higher it is, the more weight is put to the mean of the residuals and the less to the range
#'  of their corresponding x-values (see Kaishev et al., 2016, for further details).
#'
#'  The default values of \code{beta} are \code{beta = 0.5} if
#'  the response is assumed to be Gaussian, \code{beta = 0.2} if it is Poisson (or Quasipoisson),
#'   while if it is Binomial, Quasibinomial or Gamma \code{beta = 0.1}, which reflect our experience of
#'   running GeDS for different underlying functional dependencies.
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
#'  each of the three alternative stopping rules can be found in Dimitrova et al. (2017).
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
#'
#' @references
#' Kaishev, V.K., Dimitrova, D.S., Haberman, S. and Verrall, R.J. (2016).
#' Geometrically designed, variable knot regression splines.
#' \emph{Computational Statistics}, \strong{31}, 1079--1105. \cr
#' DOI: \href{https://doi.org/10.1007/s00180-015-0621-7}{doi.org/10.1007/s00180-015-0621-7}
#'
#' Dimitrova, D.S., Kaishev, V.K., Lattuada A. and Verrall, R.J. (2017).
#' Geometrically designed, variable knot splines in Generalized (Non-)Linear Models.
#' Available at \href{http://openaccess.city.ac.uk/18460/}{openaccess.city.ac.uk}
#'
#'
#' @seealso \code{\link{NGeDS}}; \code{\link{GeDS-Class}}; S3 methods such as \code{\link{coef.GeDS}}, \code{\link{deviance.GeDS}},
#'  \code{\link{knots.GeDS}}, \code{\link{print.GeDS}} and \code{\link{predict.GeDS}}; \code{\link{Integrate}}
#'   and \code{\link{Derive}}; \code{\link{PPolyRep}}.
#'
#' @examples
#' ######################################################################
#' # Generate a data sample for the response variable Y and the covariate
#' # X assuming Poisson distributed error and log link function
#' # See section 4.1 in Dimitrova et al. (2017)
#' set.seed(123)
#' N <- 500
#' f_1 <- function(x) (10*x/(1+100*x^2))*4+4
#' X <- sort(runif(N, min = -2, max = 2))
#' # Specify a model for the mean of Y to include only a component
#' # non-linear in X, defined by the function f_1
#' means <- exp(f_1(X))
#' # Generate Poisson distributed Y according to the mean model
#' Y <- rpois(N, means)
#'
#' # Fit a Poisson GeDS regression using GGeDS
#' (Gmod <- GGeDS(Y ~ f(X), beta = 0.2, phi = 0.995, family = poisson(),
#'                 Xextr = c(-2,2)))
#'
#' # Plot the quadratic and cubic GeDS fits
#' plot(X,log(Y),xlab = "x", ylab = expression(f[1](x)))
#' lines(Gmod, n = 3, col = "red")
#' lines(Gmod, n = 4, col = "blue", lty = 2)
#' legend("topleft", c("Quadratic", "Cubic"),
#'        col = c("red", "blue"), lty = c(1,2))
#'
#' # Generate GeDS prediction at X=0, first on the response scale and then on
#' # the predictor scale
#' predict(Gmod, n = 3, newdata = data.frame(X = 0))
#' predict(Gmod, n = 3, newdata = data.frame(X = 0), type = "link")
#'
#' # Apply some of the other available methods, e.g.
#' # knots, coefficients and deviance extractions for the
#' # quadratic GeDS fit
#' knots(Gmod)
#' coef(Gmod)
#' deviance(Gmod)
#'
#' # the same but for the cubic GeDS fit
#' knots(Gmod, n = 4)
#' coef(Gmod, n = 4)
#' deviance(Gmod, n = 4)
#'
#'
#' ##########################################
#' # A real data example
#' # See Dimitrova et al. (2017), Section 4.2
#'
#' data("coalMining")
#' (Gmod2 <- GGeDS(formula = accidents ~ f(years), beta = 0.1, phi = 0.98,
#'                  family = poisson(), data = coalMining))
#' (Gmod3 <- GGeDS(formula = accidents ~ f(years), beta = 0.1, phi = 0.985,
#'                  family = poisson(), data = coalMining))
#' plot(coalMining$years, coalMining$accidents, type = "h", xlab = "Years",
#'      ylab = "Accidents")
#' lines(Gmod2, tr = exp, n = 4, col = "red")
#' lines(Gmod3, tr = exp, n = 4, col = "blue", lty = 2)
#' legend("topright", c("phi = 0.98","phi = 0.985"), col = c("red", "blue"),
#'        lty=c(1, 2))
#'
#'
#' \dontrun{
#' ##########################################
#' # The same regression in the example of GeDS
#' # but assuming Gamma and Poisson responses
#' # See Dimitrova et al. (2017), Section 4.2
#'
#' data('BaFe2As2')
#' (Gmod4 <- GGeDS(intensity ~ f(angle), data = BaFe2As2, beta = 0.6, phi = 0.995, q = 3,
#'                 family = Gamma(log), stoptype = "RD"))
#' plot(Gmod4)
#'
#' (Gmod5 <- GGeDS(intensity ~ f(angle), data = BaFe2As2, beta = 0.1, phi = 0.995, q = 3,
#'                 family = poisson(), stoptype = "SR"))
#' plot(Gmod5)
#' }
#'
#' ##########################################
#' # Life tables
#' # See Dimitrova et al. (2017), Section 4.2
#'
#' data(EWmortality)
#' attach(EWmortality)
#' (M1 <- GGeDS(formula = Deaths ~ f(Age) + offset(log(Exposure)),
#'               family = poisson(), phi = 0.99, beta = 0.1, q = 3,
#'               stoptype = "LR"))
#'
#' Exposure_init <- Exposure + 0.5 * Deaths
#' Rate <- Deaths / Exposure_init
#' (M2 <- GGeDS(formula = Rate ~ f(Age), weights = Exposure_init,
#'               family = quasibinomial(), phi = 0.99, beta = 0.1,
#'               q = 3, stoptype = "LR"))
#'
#'
#' op <- par(mfrow=c(2,2))
#' plot(Age, Deaths/Exposure, ylab = expression(mu[x]), xlab = "Age")
#' lines(M1, n = 3, tr = exp, lwd = 1, col = "red")
#' plot(Age, Rate, ylab = expression(q[x]), xlab = "Age")
#' lines(M2, n = 3, tr = quasibinomial()$linkinv, lwd = 1, col = "red")
#' plot(Age, log(Deaths/Exposure), ylab = expression(log(mu[x])), xlab = "Age")
#' lines(M1, n = 3, lwd = 1, col = "red")
#' plot(Age, quasibinomial()$linkfun(Rate), ylab = expression(logit(q[x])), xlab = "Age")
#' lines(M2, n = 3, lwd = 1, col = "red")
#' par(op)
#'
GGeDS<- function(formula, data, family = gaussian(), weights, beta,  phi = 0.99,
                 min.intknots, max.intknots, q = 2L, Xextr = NULL, show.iters = FALSE, stoptype = "SR"){
  save <- match.call()
  if (missing(data))
    data <- environment(formula)
  wn <- match("weights",names(save),0L)
  wn <- save[c(1L,wn)]
  wn[[1L]] <- quote(list)
  weights <- eval(wn,data,parent.frame())$weights

  if(missing(min.intknots)) min.intknots <- 0
  newdata <- read.formula(formula, data)
  X <- newdata$X
  Y <- newdata$Y
  Z <- newdata$Z
  ncz <- if(is.null(Z)) 0 else NCOL(Z)
  if(missing(max.intknots)) max.intknots <- length(unique(X)) - 2 - ncz

  offset <- newdata$offset
  if(!is.null(Z)) Z <- as.matrix(Z)
  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
    family <- family()
  if (is.null(family$family)) {
    stop("'family' not recognized")
  }
  if(missing(beta)){
    beta <- 0.5
    if(family$family=="gaussian") beta <- 0.5
    if(family$family%in%c("poisson","quasipoisson")) beta <- 0.2
    if(family$family%in%c("binomial","quasibinomial")) beta <- 0.1
    if(family$family=="Gamma") beta <- 0.1
  }
  if(is.null(weights)) weights=rep(1,NROW(X))
  weights <- as.numeric(weights)
  if(is.null(Z)) {
    if(NROW(X)!=length(Y) || NROW(X)!=length(weights)) stop("length of 'X', 'Y' and 'weights' must match")
  } else {
    if(NROW(X)!=length(Y) || NROW(X)!=NROW(Z) || NROW(X)!=length(weights)) stop("length of 'X', 'Y', 'Z' and 'weights' must match")
  }
  if(is.null(X) || is.null(Y)) stop("Null arguments cannot be passed to this function")
  if(length(phi) != 1 || length(beta) != 1  || length(max.intknots) != 1 ||
     length(q) != 1 || length(show.iters) != 1)
    stop("'phi', 'beta', 'max.intknots', 'q' and 'show.iters' must have length 1")
  if(beta > 1 || beta < 0) stop("'beta' should be a value in [0,1]")
    phi <- as.numeric(phi)
    if(phi > 1 || phi < 0) stop("'phi' should be a value in [0,1]")


  if(ncol(X)==1){
    max.intknots <- as.integer(max.intknots)
    min.intknots <- as.integer(min.intknots)
    q <- as.integer(q)
    show.iters <- as.logical(show.iters)
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
        Z <- Z[!tmp,]
        offset <- offset[!tmp]
        weights <- weights[!tmp]
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
    if(length(Xextr) !=2) stop("'extr' must have length 2")
    out <- GenUnivariateFitter(
      X=X, Y=Y, Z=Z, offset=offset,
      weights = weights, family=family,
      beta=beta, phi = phi, min.intknots = min.intknots,
      max.intknots = max.intknots,
      q = q, extr=Xextr, show.iters=show.iters,stoptype=stoptype)
  } else {
    if(ncol(X)==2){
      stop("Bivariate version not yet implemented")
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
