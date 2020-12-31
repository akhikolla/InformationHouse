read.formula <- function(formula,data,weights,offset){
  mt <- if (missing(data)){
    terms(formula, "f")
  } else terms(formula, "f", data = data)
  if(attr(mt,"intercept")==1) {
    attr(mt,"intercept") <- 0
    #warning("Intercept will be included in the basis functions")
  } else {
    warning("An intercept will be included in the basis functions")
  }
#  data <- new.env(parent=as.environment(data))
  data <- as.list(data)
  data$f <- function(x,xx=NULL,...) {
    if(!missing(...)) stop("Algorithm supports at most two variables in 'f'")
    cbind(x,xx)
  }
  spec <- attr(mt,"specials")$f

  if(length(spec)!=1) stop("Formula incorrectly specified. Read documentation for further informations")
  mm <- model.matrix(mt,data)
  mf <- model.frame(mt,data,na.action=NULL)
  Y <- model.response(mf, type="numeric")
  attr(Y,"names")<- NULL
  #names(mf) <- NULL
  X <- mf[,spec]
  if(ncol(mm)>ncol(X)) {
    Z <- mf[,-c(spec,attr(mt,"response"),attr(mt,"offset")),drop=T]
  } else {
    Z <- NULL
  }

  offset <- rep(0, nrow(X))
  if (!is.null(off.num <- attr(mt, "offset")))
    for (i in off.num) offset <- offset + eval(attr(mt,"variables")[[i + 1]], data)



  out<-list("X" = X, "Y" = Y, "Z" = Z, "offset" = offset, "terms" = mt, "model.matrix" = mm)
  return(out)
}



#' Formula for the predictor model
#'
#' @description
#' A description of the structure of a predictor model to be fitted using  \code{\link{NGeDS}}
#' and/or \code{\link{GGeDS}} and how this information can be extracted from a \code{\link{GeDS-class}} object.
#'
#' @aliases formula.GeDS
#' @export
#'
#' @param x Fitted \code{\link{GeDS-class}} object, tipically produced by
#' \code{\link{NGeDS}} or \code{\link{GGeDS}} from which the predictor model \code{\link[stats]{formula}} should be extracted.
#' @param ... Unused in this case.
#'
#' @details
#' In the GeDS GNM (GLM) regression, implemented in \code{\link{NGeDS}}
#' and \code{\link{GGeDS}}, it is assumed that the mean of the response variable transformed using an
#' appropriate link function is modelled through a possibly multivariate predictor model involving two components:
#' a GeD variable knot spline regression component involving up to two of the
#'  independent variables and a parametric component with respect to  the remaining independent variables.
#'  The formula is used to specify the structure of such a possibly multivariate predictor model.
#'
#'  The formulae that are input in \code{\link{NGeDS}} and \code{\link{GGeDS}} are similar to those input
#'  in \code{\link[stats]{lm}} or \code{\link[stats]{glm}} except that the function \code{\link{f}} should be
#'  specified in order to identify which of the covariates enter the GeD spline regression part
#'  of the predictor model. For example,  if the predictor model is univariate and it links the transformed means of \code{y}
#'  to \code{x1}, the predictor has only a GeD spline component and the \code{\link[stats]{formula}}
#'  should be in the form \code{y ~ f(x1)}.
#'
#'  As noted, there may be additional independent variables, \code{x2}, \code{x3}, ... which may
#' enter linearly into the parametric component of the predictor model and not be part of the
#' GeD spline regression component. For example one may use the formula
#' \code{y ~ f(x1) + x2 + x3} which assumes a spline regression only between the transformed mean of \code{y}
#' and \code{x1}, while \code{x2} and \code{x3} enter the predictor model just linearly.
#'
#' In the current version of the package, \code{\link{GGeDS}} is univariate, therefore only one covariate
#' which enters the spline regression component can be specified.
#'
#' In contrast, the function \code{\link{NGeDS}}, generates also bivariate GeDS regression models.
#' Therefore, if the functional dependence of the mean of the response variable \code{y} on \code{x1} and
#'  \code{x2} needs to be jointly modelled and there are no other covariates, the formula for the corresponding
#'  two dimensional predictor model should be specified as \code{y ~ f(x1,x2)}.
#'
#'  Within the argument \code{formula}, similarly as in other \R functions, it is possible to
#'  specify one or more offset variables, i.e. known terms with fixed regression coefficients equal to 1.
#'  These terms should be identified via the function \code{\link[stats]{offset}}.
#'
formula.GeDS <- function(x,...){
  frm <- x$Formula
  if(is.null(frm)) stop("Unable to extract the formula. \n
                        Please re-fit using 'NGeDS' or 'GGeDS'")
  frm
}


# #' versions of GeDS have been implemented correspondingly
# #' in \code{\link{GGeDS}} and \code{\link{NGeDS}}.
# #' If more than two covariates are
# #' specifying these arguments will return an error.



#' Defining the covariates for the spline component in a GeDS formula.
#'
#' @description
#' In general the GeDS predictor model may include a GeD spline regression component
#' with respect to part of the independent variables and a parametric component
#' in which the remaining covariates may enter as additive terms.
#'
#' The function \code{f} is to be used in the \code{\link[=formula.GeDS]{formula}} argument of \code{\link{NGeDS}}
#' or \code{\link{GGeDS}} in order to specify which independent variables (covariates)
#' should be included in the GeD spline regression component of the predictor model.
#'
#' @note This function is intended to be used only as part of the \code{\link[=formula.GeDS]{formula}}
#' in a GeDS regression via \code{\link{NGeDS}} or \code{\link{GGeDS}}
#' and not to be called in other cases by the user.
#'
#' @param x numeric vector containing \eqn{N} sample values of the covariate chosen to enter the spline
#' regression component of the predictor model.
#' @param xx numeric vector containing \eqn{N} sample values for the second covariate
#' (in case \code{\link{NGeDS}} is run for two dimensions).
#' It has to be either \code{NULL} (the default) or a vector of size \eqn{N}, same as \code{x}.
#' @param ... further arguments. As GeDS currently allows for up to two covariates,
#' specification of further arguments will return an error.
#'
#' @examples
#' # Generate a data sample for the response variable Y and
#' # the covariates X, reg1, reg2 and off
#' set.seed(123)
#' N <- 500
#' f_1 <- function(x) (10*x/(1+100*x^2))*4+4
#' X <- sort(runif(N ,min = -2, max = 2))
#' reg1 <- runif(500, min = -0.1, max = 0.1)
#' reg2 <- runif(500, min = -0.2, max = 0.2)
#' off <- runif(500, min = -1, max = 1)
#' # Specify a model for the mean of Y to include a component non linear
#' # in X defined by the function f_1 and a linear one in the other covariates
#' means <- f_1(X) + 2*reg1 + 0.5*reg2 + off
#' # Add Normal noise to the mean of Y
#' Y <- rnorm(N, means, sd = 0.1)
#'
#' # Specify a formula that will be used to model Y as a
#' # function of X, reg1, reg2 and off.
#' # The covariate X is for the spline component modeled as GeDS,
#' # reg1 and reg2 enter linearly, off is an offset, i.e. no coefficient
#' # will be estimated for it
#' formula <- Y ~ f(X) + reg1 + reg2 + offset(off)
#'
#' # Fit a GeDS model specified in formula using NGeDS
#' (Gmod <- NGeDS(formula, beta = 0.6, phi = 0.995, Xextr = c(-2,2)))
#'
#'
#' @seealso \link{NGeDS}; \link{GGeDS}.
f <- function(x,xx=NULL,...) {
  if(!missing(...)) stop("Algorithm supports at most two variables in 'f'")
  cbind(x,xx)
}


# this is to get the names of the Z variable(s)
getZnames <- function(out){
  names <- colnames(out$model.matrix)
  id <- attr(out$model.matrix,"assign")
  spec <- attr(out$terms,"specials")$f-1
  znames <- names[id!=spec]
  znames
}


