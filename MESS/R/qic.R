#' Quasi Information Criterion
#'
#' Function for calculating the quasi-likelihood under the independence model
#' information criterion (QIC), quasi-likelihood, correlation information
#' criterion (CIC), and corrected QIC for one or several fitted geeglm model
#' object from the geepack package.
#'
#' QIC is used to select a correlation structure. The QICu is used to compare
#' models that have the same working correlation matrix and the same
#' quasi-likelihood form but different mean specifications. CIC has been
#' suggested as a more robust alternative to QIC when the model for the mean
#' may not fit the data very well and when models with different correlation
#' structures are compared.
#'
#' Models with smaller values of QIC, CIC, QICu, or QICC are preferred.
#'
#' If the MASS package is loaded then the \code{\link{ginv}} function is used
#' for matrix inversion. Otherwise the standard \code{\link{solve}} function is
#' used.
#'
#' @aliases QIC QIC.geeglm QIC.geekin QIC.ordgee
#' @param object a fitted GEE model from the geepack package. Currently only
#' works on geeglm objects
#' @param tol the tolerance used for matrix inversion
#' @param \dots optionally more fitted geeglm model objects
#' @return A vector or matrix with the QIC, QICu, quasi likelihood, CIC, the
#' number of mean effect parameters, and the corrected QIC for each GEE object
#' @author Claus Ekstrom \email{claus@@rprimer.dk}
#' @seealso \code{geeglm}
#' @references Pan, W. (2001). \emph{Akaike's information criterion in
#' generalized estimating equations}. Biometrics, 57, 120-125.\cr Hardin, J.W.
#' and Hilbe, J.M. (2012). \emph{Generalized Estimating Equations, 2nd
#' Edition}, Chapman and Hall/CRC: New York. \cr Hin, L.-Y. and Wang, Y-G.
#' (2009). \emph{Working-correlation-structure identification in generalized
#' estimating equations}, Statistics in Medicine 28: 642-658. \cr Thall, P.F.
#' and Vail, S.C. (1990). \emph{Some Covariance Models for Longitudinal Count
#' Data with Overdispersion}.  Biometrics, 46, 657-671.
#' @keywords htest
#' @examples
#'
#' library(geepack)
#' data(ohio)
#' fit <- geeglm(resp ~ age + smoke + age:smoke, id=id, data=ohio,
#'              family=binomial, corstr="exch", scale.fix=TRUE)
#' QIC(fit)
#'
#' @rdname QIC
#' @export
QIC.geeglm <- function(object, tol=.Machine$double.eps, ...) {

  #
  # The majority of this code was taken from the internet
  # I added a bit of functionality and made the whole interface smoother

  if (! ("geeglm" %in% class(object)) ) {
    stop("QIC requires a geeglm object as input")
  }

  # Setup functions
  invert <- if ("MASS" %in% loadedNamespaces()) {
    MASS::ginv
  } else { solve }

  # Missing:
  # Check correct handling of link and family functions

  # Create function to make the computations
  computeqic <- function(object) {
    # Fitted and observed values for quasi likelihood
    mu <- object$fitted.values
    y  <- object$y

    # Quasi Likelihood for Poisson
    # quasi.R <- sum((y*log(mu.R)) - mu.R) # poisson()$dev.resids - scale and weights = 1
    type <- family(object)$family
    quasi <- switch(type,
                    poisson = sum((y*log(mu)) - mu),
                    gaussian = sum(((y - mu)^2)/-2),
                    binomial = sum(y*log(mu/(1 - mu)) + log(1 - mu)),
                    Gamma = sum(-y/(mu - log(mu))),
                    stop("Error: distribution not recognized"))

    # Fit model with independence correlation structure
    object$call$corstr <- "independence"
    object$call$zcor <- NULL
    model.indep <- eval(object, parent.frame())
    # model.indep <- update(object, corstr="independence",zcorr=NULL)

    # Trace term (penalty for model complexity)
    AIinverse <- invert(model.indep$geese$vbeta.naiv, tol=tol)
    Vr <- object$geese$vbeta
    trace <- sum(diag(AIinverse %*% Vr))
    params <- length(coef(object)) # Mean parameters in the model

    kpm <- params+length(object$geese$alpha)

    # QIC
    QIC <- -2*(quasi - trace)
    QICu <- -2*(quasi - params)
    QICC <- QIC + (2*kpm*(kpm+1))/(length(object$residuals)-kpm-1)
    output <- c(QIC, QICu, quasi, trace, params, QICC)
    names(output) <- c("QIC", "QICu", "Quasi Lik", "CIC", "params", "QICC")
    output
  }

  if (length(list(...))) {
    # Make the computations
    results <- lapply(list(object, ...), computeqic)

    # Check same data size
    check <- sapply(list(object, ...), function(x) {
      length(x$y)
    })

    if (any(check != check[1]))
      warning("models are not all fitted to the same number of observations")

    # Merge the results together in a data.matrix
    res <- do.call("rbind", results)

    # Set the row names corresponding to the models
    Call <- match.call()
    Call$k <- NULL
    row.names(res) <- as.character(Call[-1L])
    res
  } else {
    computeqic(object)
  }
}



#' @rdname QIC
#' @export
QIC.ordgee <- function(object, tol = .Machine$double.eps, ...) {

  #
  # The majority of this code was taken from the internet
  # I added a bit of functionality and made the whole interface smoother

  if (! ("geeglm" %in% class(object)) ) {
    stop("QIC requires a geeglm object as input")
  }

  # Setup functions
  invert <- if ("MASS" %in% loadedNamespaces()) {
    MASS::ginv
  } else { solve }

  # Missing:
  # Check correct handling of link and family functions

  # Create function to make the computations
  computeqic <- function(object) {
    # Fitted and observed values for quasi likelihood
    mu <- object$fitted.values
    y  <- object$y

    # Quasi Likelihood for Poisson
    # quasi.R <- sum((y*log(mu.R)) - mu.R) # poisson()$dev.resids - scale and weights = 1
    type <- family(object)$family
    quasi <- switch(type,
                    poisson = sum((y*log(mu)) - mu),
                    gaussian = sum(((y - mu)^2)/-2),
                    binomial = sum(y*log(mu/(1 - mu)) + log(1 - mu)),
                    Gamma = sum(-y/(mu - log(mu))),
                    stop("Error: distribution not recognized"))

    # Fit model with independence correlation structure
    object$call$corstr <- "independence"
    object$call$zcor <- NULL
    model.indep <- eval(object, parent.frame())
    # model.indep <- update(object, corstr="independence",zcorr=NULL)

    # Trace term (penalty for model complexity)
    AIinverse <- invert(model.indep$geese$vbeta.naiv)
    Vr <- object$geese$vbeta
    trace <- sum(diag(AIinverse %*% Vr))
    params <- length(coef(object)) # Mean parameters in the model

    kpm <- params+length(object$geese$alpha)

    # QIC
    QIC <- -2*(quasi - trace)
    QICu <- -2*(quasi - params)
    QICC <- QIC + (2*kpm*(kpm+1))/(length(object$residuals)-kpm-1)
    output <- c(QIC, QICu, quasi, trace, params, QICC)
    names(output) <- c("QIC", "QICu", "Quasi Lik", "CIC", "params", "QICC")
    output
  }

  if (length(list(...))) {
    # Make the computations
    results <- lapply(list(object, ...), computeqic)

    # Check same data size
    check <- sapply(list(object, ...), function(x) {
      length(x$y)
    })

    if (any(check != check[1]))
      warning("models are not all fitted to the same number of observations")

    # Merge the results together in a data.matrix
    res <- do.call("rbind", results)

    # Set the row names corresponding to the models
    Call <- match.call()
    Call$k <- NULL
    row.names(res) <- as.character(Call[-1L])
    res
  } else {
    computeqic(object)
  }
}



## QIC.binomial <- function(object, ...) {

##   #
##   # The majority of this code was taken from the internet
##   # I added a bit of functionality and made the whole interface smoother

##   if (! ("geese" %in% class(object)) ) {
##     stop("QIC requires a geese object as input")
##   }

##   # Setup functions
##   invert <- if ("MASS" %in% loadedNamespaces()) {
##     MASS:::ginv
##   } else { solve }

##   # Missing:
##   # Check correct handling of link and family functions

##   # Create function to make the computations
##   computeqic <- function(object) {
##     # Fitted and observed values for quasi likelihood
##     # Compute the linear predictor
##     glmcall <- object$call
##     glmcall$id <- glmcall$jack <- glmcall$control <- glmcall$corstr <- glmcall$waves <- glmcall$zcor <- glmcall$std.err <- glmcall$scale.fix <- glmcall$scale.value<- glmcall$z <- glmcall$family <- NULL
##     glmcall[[1]] <- as.name("model.frame")
##     mf <- eval(glmcall, parent.frame())
##     X <- model.matrix(formula(object), data=mf)

##     N <- nrow(X)
##     offset <- model.offset(mf)
##     if (is.null(offset))
##         offset <- rep(0, N)

##     if (is.null(object$call$offset))
##         mu <- as.vector(X %*% object$beta) else mu <- offset + X %*% object$beta

##     y  <- as.numeric(levels(mf[,1]))[mf[,1]]
##     if (length(unique(y))!=2)
##         stop("QIC.binomial only works for binary data")

##     # Quasi Likelihood for Poisson
##     # quasi.R <- sum((y*log(mu.R)) - mu.R) # poisson()$dev.resids - scale and weights = 1
## #    type <- family(object)$family
##     type <- "binomial"
##     quasi <- switch(type,
##                     poisson = sum((y*log(mu)) - mu),
##                     gaussian = sum(((y - mu)^2)/-2),
##                     binomial = sum(y*mu + log(1 - exp(mu)/(1+exp(mu)))),
##                     Gamma = sum(-y/(mu - log(mu))),
##                     stop("Error: distribution not recognized"))

##     # Fit model with independence correlation structure
##     object$call$corstr <- "independence"
##     object$call$zcor <- object$call$z <- NULL
##     model.indep <- eval(object$call, parent.frame())

##     # model.indep <- update(object, corstr="independence")

##     # Trace term (penalty for model complexity)
##     AIinverse <- invert(model.indep$vbeta.naiv)
##     Vr <- object$vbeta
##     trace <- sum(diag(AIinverse %*% Vr))
##     params <- length(object$beta) # Mean parameters in the model

##     kpm <- params+length(object$alpha)

##     # QIC
##     QIC <- -2*(quasi - trace)
##     QICu <- -2*(quasi - params)
##     QICC <- QIC - (2*kpm*(kpm+1))/(length(object$residuals)-kpm-1)
##     output <- c(QIC, QICu, quasi, trace, params, QICC)
##     names(output) <- c("QIC", "QICu", "Quasi Lik", "CIC", "params", "QICC")
##     output
##   }

##   if (length(list(...))) {
##     # Make the computations
##     results <- lapply(list(object, ...), computeqic)

##     # Check same data size
##     check <- sapply(list(object, ...), function(x) {
##       n <- length(x$y)
##     })

##     if (any(check != check[1]))
##       warning("models are not all fitted to the same number of observations")

##     # Merge the results together in a data.matrix
##     res <- do.call("rbind", results)

##     # Set the row names corresponding to the models
##     Call <- match.call()
##     Call$k <- NULL
##     row.names(res) <- as.character(Call[-1L])
##     res
##   } else {
##     computeqic(object)
##   }
## }


#' @rdname QIC
#' @export
QIC.geekin <- function(object,  tol = .Machine$double.eps, ...) {

  # This functions is only needed to replace class
  # geeglm to make sure the regular
  # QIC function works

  if (! ("geeglm" %in% class(object)) ) {
    stop("QIC requires a geekin object as input")
  }

  object$call[[1]] <- as.name("geeglm")
  object$call$varlist <- NULL
  object$call$na.action <- NULL

  # Swap class around
  class(object) <- c("geeglm", unique(class(object)))
  QIC(object)
}


#' @rdname QIC
#' @export
QIC <- function(object,  tol = .Machine$double.eps, ...) {
  UseMethod("QIC")
}
