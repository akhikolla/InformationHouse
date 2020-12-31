#' Coef method for GeDS objects
#' @export
#' @rdname coef
#' @aliases coef.GeDS coefficients.GeDS
#'
#' @description Methods for the functions \code{\link[stats]{coef}} and \code{\link[stats]{coefficients}}
#'  that allow to extract the
#' estimated coefficients of a fitted GeDS regression from a \code{\link{GeDS-Class}} object.
#' @param object the  \code{\link{GeDS-class}} object from which the
#' coefficients of the selected GeDS regression should be extracted.
#' @param n integer value (2, 3 or 4) specifying the order (\eqn{=} degree \eqn{+ 1}) of the GeDS fit
#' whose coefficients should be extracted.
#' By default equal to \code{3L}.
#' Non-integer values will be passed to the function \code{\link{as.integer}}.
#' @param onlySpline logical variable specifying whether only the coefficients for the GeDS  component
#' of the fitted multivariate regression model should be extracted or alternatively also the coefficients
#'  of the parametric component should also be extracted.
#' @param ... potentially further arguments (required by the definition of the generic function).
#' They will be ignored, but with a warning.
#'
#' @details
#' These are simple methods for the functions \code{\link[stats]{coef}} and \code{\link[stats]{coefficients}}.
#'
#' As \code{\link{GeDS-class}} objects contain three different fits (linear, quadratic and cubic), it is possible
#' to specify the order of the fit for which GeDS regression coefficients are required via the input argument \code{n}.
#'
#' As mentioned in the details of \code{\link[=formula.GeDS]{formula}}, the predictor model may be multivariate and it
#' may include a GeD
#' spline component whereas the remaining variables may be part of a
#' parametric component. If the \code{onlySpline} argument is set to \code{TRUE} (the default value),
#' only the coefficients corresponding to the GeD spline component of order \code{n} of the multivariate
#' predictor model are extracted.
#'
#'
#' @return
#' A named vector containing the required coefficients of the fitted multivariate predictor model.
#' The coefficients corresponding to the variables that enter the
#'  parametric component of the fitted multivariate predictor model
#' are named as the variables themselves. The  coefficients of the GeDS component
#' are coded as "\code{N}" followed by the index of the corresponding B-spline.
#'
#' @seealso \code{\link[stats]{coef}} for the standard definition; \code{\link{NGeDS}} for examples.
#' @examples
#' # Generate a data sample for the response variable
#' # and the covariates
#' set.seed(123)
#' N <- 500
#' f_1 <- function(x) (10*x/(1+100*x^2))*4+4
#' X <- sort(runif(N ,min = -2, max = 2))
#' Z <- runif(N)
#' # Specify a model for the mean of the response Y to be a superposition of
#' # a non-linear component f_1(X), a linear component 2*Z and a
#' # free term 1, i.e.
#' means <- f_1(X) + 2*Z + 1
#' # Add normal noise to the mean of y
#' Y <- rnorm(N, means, sd = 0.1)
#'
#' # Fit to this sample a predictor model of the form f(X) + Z, where
#' # f(X) is the GeDS component and Z is the linear (additive) component
#' # see ?formula.GeDS for details
#' (Gmod <- NGeDS(Y ~ f(X) + Z, beta = 0.6, phi = 0.995, Xextr = c(-2,2)))
#'
#' # Extract the GeD spline regression coefficients
#' coef(Gmod, n = 3)
#'
#' # Extract all the coefficients, including the one for the linear component
#' coef(Gmod, onlySpline = FALSE, n = 3)
coef.GeDS <- function(object, n = 3L, onlySpline = TRUE, ...){
  if(!missing(...)) warning("Only 'object', 'n' and 'onlySpline' arguments will be considered")
  n <- as.integer(n)
  if(!(n %in% 2L:4L)) {
    n <- 3L
    warning("'n' incorrectly specified. Set to 3.")
  }
  if(n == 3L){
    theta <- object$Quadratic$Theta
    nth <- length(object$Quadratic$Polygon$Kn)
  }
  if(n == 4L){
    theta <- object$Cubic$Theta
    nth <- length(object$Cubic$Polygon$Kn)
  }
  if(n == 2L){
    theta <- object$Linear$Theta
    nth <- length(object$Linear$Polygon$Kn)
  }
  if(!is.null(object$Args$Z) & !onlySpline){
    znames <- attr(object$terms,"term.labels")[-1]
    names(theta) <- c(paste0("N",1:nth),znames)
  } else {
    theta <- theta[1:nth]
    names(theta) <- paste0("N",1:nth)
  }
  return(theta)
}
#' @rdname coef
#' @export
#'
coefficients.GeDS <- function(object, n = 3L, onlySpline = TRUE, ...){
  coef.GeDS(object = object, n=n, onlySpline = onlySpline, ...)
}

#' Deviance method for GeDS objects
#' @export
#' @aliases deviance.GeDS
#' @description Method for the function \code{\link[stats]{deviance}} that allows the user to extract  the value
#' of the deviance corresponding to a selected GeDS fit from a \code{\link{GeDS-Class}} object.
#' @param object the \code{\link{GeDS-class}} object from which the deviance should be extracted.
#' @param n integer value (2, 3 or 4) specifying the order (\eqn{=} degree \eqn{+ 1}) of the GeDS fit
#' whose deviance should be extracted.
#' By default equal to \code{3L}.
#' Non-integer values will be passed to the function \code{\link{as.integer}}.
#' @param ... potentially further arguments (required by the definition of the generic function).
#' They will be ignored, but with a warning.
#'
#' @details This is a method for the function \code{\link[stats]{deviance}}.
#' As \code{\link{GeDS-class}} objects contain three different fits (linear, quadratic and cubic),
#' it is possible
#' to specify the order of the GeDS fit for which the deviance is required via the input argument \code{n}.
#'
#' @seealso \code{\link[stats]{deviance}} for the standard definition; \code{\link{GGeDS}} for examples.
#' @return A numeric value corresponding to the  deviance of the selected GeDS fit.
deviance.GeDS <- function(object, n = 3L, ...){
  if(!missing(...)) warning("Only 'object' and 'n' arguments will be considered")
  if(length(n)!=1) stop("Only one Deviance at each time can be extracted")
  n <- as.integer(n)
  if(!(n %in% 2L:4L)) {
    n <- 3L
    warning("'n' incorrectly specified. Set to 3.")
  }
  if(n == 3L){
    dev <- object$Quadratic$RSS
  }
  if(n == 4L){
    dev <- object$Cubic$RSS
  }
  if(n == 2L){
    dev <- object$Linear$RSS
  }
  dev <- as.numeric(dev)
  return(dev)
}



#' Predict method for GeDS objects
#' @aliases predict.GeDS
#' @export
#' @description  This is a user friendly method to compute predictions from GeDS objects.
#' @param object the  \code{\link{GeDS-class}} object for which the
#'  computation of the predicted values is required.
#' @param newdata an optional \code{data.frame}, \code{list} or \code{environment} containing values
#'  of the independent variables for  which predicted values of the predictor model
#'  (including the GeDS and the parametric components) should be computed.
#' If left empty the values are extracted from the object \code{x} itself.
#' @param type character string indicating the type of prediction required.
#' By default it is equal to \code{"response"}, i.e. the result is on the scale of the response variable.
#' See details for the other options.
#' @param n integer value (2, 3 or 4) specifying the order (\eqn{=} degree \eqn{+ 1}) of the GeDS fit
#' whose predicted values should be computed.
#' By default equal to \code{3L}.
#' Non-integer values will be passed to the function \code{\link{as.integer}}.
#' @param ... potentially further arguments (required by the definition of the generic function.).
#' They are ignored, but with a warning.
#'
#' @details This is a method for the function \code{\link[stats]{predict}} that
#' allows the user to handle \code{\link{GeDS-Class}} objects.
#'
#' In analogy with the function \code{\link[stats]{predict.glm}} in the \pkg{stats} package,
#' the user can specify the scale on
#' which the predictions should be computed through the argument \code{type}.
#' If the predictions are required to be on the scale of the response variable, the user should set
#' \code{type = "response"}, which is the default. Alternatively
#' if one wants the predictions to be on the predictor scale, it is necessary to set \code{type = "link"}.
#'
#' By specifying \code{type = "terms"}, it is possible to inspect the predicted values
#'  separately for each single independent variable which enter either the GeD spline
#'  component or the parametric component of the predictor model.
#'  In this case the returned result is a matrix whose columns correspond to the terms supplied
#'  via \code{newdata} or extracted from the \code{object}.
#'
#' As GeDS objects contain three different fits (linear, quadratic and cubic), it is possible
#' to specify the order
#' for which GeDS predictions are required via the input argument \code{n}.
#'
#' @seealso \code{\link[stats]{predict}} for the standard definition; \code{\link{GGeDS}} for examples.
#' @return A numeric vector corresponding to the predicted values (if \code{type = "link"} or
#' \code{type = "response"}). If \code{type = "terms"} a numeric matrix with a column per term.
predict.GeDS <- function(object, newdata, type = c("response", "link", "terms"), n=3L, ...){
  if(!missing(...)) warning("Only 'object', 'newdata, 'type' and 'n' arguments will be considered")
  n <- as.integer(n)
  mt <- object$terms
  if (!inherits(object, "GeDS"))
    warning("calling predict.GeDS(<fake-GeDS-object>) ...")
  if (missing(newdata) || is.null(newdata)) {
    X <- object$Args$X
    Z <- object$Args$Z
    offset <- object$Args$offset
    if(is.null(offset)) offset <- rep(0, NROW(X))
  }
  else {
    mt <- delete.response(mt)
    newdata <- as.list(newdata)
    newdata$f <- f
    mm <- model.matrix(mt,newdata)
    mf <- model.frame(mt,newdata)
    spec <- attr(mt,"specials")$f
    X <- mf[,spec]
    if(ncol(mm)>ncol(X)) {
      Z <- mf[,-c(spec,attr(mt,"response")),drop=T]
    } else {
      Z <- NULL
    }

    offset <- rep(0, NROW(X))
    if (!is.null(off.num <- attr(mt, "offset")))
      for (i in off.num) offset <- offset +
      eval(attr(mt, "variables")[[i + 1]], newdata)
    if (!is.null(object$call$offset))
      offset <- offset + eval(object$call$offset, newdata)
  }
  kn <- knots(object,n=n,options="all")

  if(min(X)<min(kn) | max(X)>max(kn)) warning("Input values out of the boundary knots")

  matrice <- splineDesign(
    knots=kn,derivs=rep(0,length(X)),x=X,ord=n,outer.ok = T)
  type <- match.arg(type)
  if(type!="terms"){
    coefs <- coef(object,n=n, onlySpline = FALSE)
    matrice2 <- cbind(matrice,Z)
    predicted <- matrice2%*%coefs+offset

    if(type=="response" & !is.null(object$Args$family)) {
      predicted <- object$Args$family$linkinv(predicted)
    }
  } else {
    coefs <- coef(object,n=n, onlySpline = TRUE)
    coefs1 <- coef(object,n=n, onlySpline = FALSE)
    predicted <- matrice%*%coefs
    colnames(predicted) <- "Spline"
    predicted1 <- if(!is.null(Z)){
      Z*matrix(coefs1[-c(1:length(coefs))],
               ncol=length(coefs1[-c(1:length(coefs))]),nrow=NROW(Z))
      } else NULL
    if(!is.null(predicted1)) {
      predicted1 <- as.matrix(predicted1)
        colnames(predicted1) <- object$znames
    }
    predicted <- cbind(predicted,predicted1)
  }
  predicted
}





#' Print method for GeDS objects
#' @export
#' @description Method for the generic function \code{\link[base]{print}} that allows to
#'  print on screen the main information related to the fitted predictor model that can be extracted
#'   from a \code{\link{GeDS-class}} object.
#' @param x the \code{\link{GeDS-class}} object for which the main information should be printed on screen.
#' @param digits number of digits to be printed.
#' @param ... potentially further arguments (required by the definition of the generic function).
#'
#' @details This method allows to print on screen basic information related to the fitted predictor model such as the
#' function \code{call}, the number of internal knots for the linear GeDS fit and the deviances
#' for the three (linear, quadratic and cubic) fitted predictor models embedded in the \code{\link{GeDS-class}} object.
#'
#' @seealso \code{\link[base]{print}} for the standard definition.
#' @return This function returns (invisibly) the same input object, but adding the slot \code{Print}
#' that contains the three sub-slots:
#' \item{Nknots}{ the number of internal knots of the linear GeDS fit}
#' \item{Deviances}{ the deviances of the three (linear, quadratic and cubic) GeDS fits}
#' \item{Call}{ the \code{call} to the function that produced the \code{x} object}
#'
print.GeDS <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
  cat("\nCall:\n", paste(deparse(x$extcall), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  kn <- knots(x,n=2,options="int")
  names(kn) <- NULL
  RMSs <- numeric(3)
  names(RMSs) <- c("Order 2","Order 3","Order 4")
  if(is.list(kn)){
    if(length(kn[[1]])||length(kn[[2]])){
      cat(paste0("Number of internal knots of the second order spline in the X direction: ", length(kn[[1]])))
      cat("\n")
      cat(paste0("Number of internal knots of the second order spline in the Y direction: ", length(kn[[2]])))
      cat("\n")
      cat("\n")
      RMSs[1] <- x$Dev.Linear
      RMSs[2] <- if(!is.null(x$Dev.Quadratic)) x$Dev.Quadratic else NA
      RMSs[3] <- if(!is.null(x$Dev.Cubic)) x$Dev.Cubic else NA
      cat("Deviances:\n")
      print.default(format(RMSs, digits = digits), print.gap = 2L,
                    quote = FALSE)

    } else {cat("No internal knots found\n")}
    cat("\n")
    prnt <- list("Nknots" = c(length(kn[[1]]),length(kn[[2]])), "Deviances" = RMSs, "Call" = x$extcall)

  } else {

  if (length(kn)) {
    cat(paste0("Number of internal knots of the second order spline: ", length(kn)))
    cat("\n")
    RMSs[1] <- x$Dev.Linear
    RMSs[2] <- if(!is.null(x$Dev.Quadratic)) x$Dev.Quadratic else NA
    RMSs[3] <- if(!is.null(x$Dev.Cubic)) x$Dev.Cubic else NA
    cat("Deviances:\n")
    print.default(format(RMSs, digits = digits), print.gap = 2L,
                  quote = FALSE)
  }
  else cat("No internal knots found\n")
  cat("\n")
  prnt <- list("Nknots" = length(kn), "Deviances" = RMSs, "Call" = x$extcall)
  }
  x$Print <- prnt
  invisible(x)
}

#' Knots method for GeDS objects
#' @export
#' @rdname knots
#' @aliases knots.GeDS
#'
#' @description Method for the generic function \code{\link[stats]{knots}} that allows the user
#' to extract vector of knots of a GeDS fit of a specified order
#' contained in a \code{\link{GeDS-class}} object.
#'
#' @param Fn the \code{\link{GeDS-class}} object from which the vector of knots for the specified GeDS fit
#' should be extracted.
#' @param n integer value (2, 3 or 4) specifying the order (\eqn{=} degree \eqn{+ 1}) of the GeDS fit
#' whose knots should be extracted.
#' By default equal to \code{3L}.
#' Non-integer values will be passed to the function \code{\link{as.integer}}.
#' @param options a character string specifying whether "\code{all}" knots, including
#' the left-most and the right-most limits of the interval embedding the observations (the default) or
#' only the "\code{internal}" knots should be extracted.
#' @param ... potentially further arguments (required for compatibility with the definition of
#' the generic function). Currently ignored, but with a warning.
#' @details This is a method for the function \code{\link[stats]{knots}} in the \pkg{stats} package.
#'
#' As \code{\link{GeDS-class}} objects contain three different fits (linear, quadratic and cubic), it is possible
#' to specify the order of the GeDS fit  whose knots  are required via the input argument \code{n}.
#'
#'
#' @seealso \code{\link[stats]{knots}} for the definition of the generic function; \code{\link{NGeDS}} and \code{\link{GGeDS}} for examples.
#'
#' @return A vector in which each element represents a knot of the GeDS fit of the required order.
#'
#'
knots.GeDS <-  function(Fn, n = 3L, options = c("all","internal"), ...){
  if(!missing(...)) warning("Arguments other than 'Fn', 'n' and 'options' currenly igored. \n Please check if the input parameters have been correctly specified.")
  options <- match.arg(options)
  n <- as.integer(n)
  if(!(n %in% 2L:4L)) {
    n <- 3L
    warning("'n' incorrectly specified. Set to 3.")
  }
  if(n == 3L){
    kn <- Fn$Quadratic.Knots
  }
  if(n == 4L){
    kn <- Fn$Cubic.Knots
  }
  if(n == 2L){
    kn <- Fn$Linear.Knots
  }
  extr <- if(options=="all") Fn$Args$extr else NULL
  if(Fn$Type!="LM - biv"){
    kn <- sort(c(rep(extr,n),kn))
  }
  return(kn)
}

# splineKnots.GeDS <-  function(object, n = 3, options = c("all","internal")){
#   knots.GeDS(Fn=object, n = n, options = options)
# }
# would be nice, but unfortunately the generic definition in pkg splines
# does not allow for ... argument. MAybe it will change in future R releases
