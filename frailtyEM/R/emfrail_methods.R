#' @export
#' @keywords internal
coef.emfrail <- function(object, ...) {
  object$coefficients
}

#' @export
#' @keywords internal
vcov.emfrail <- function(object, type = c("regular", "adjusted"), ...) {

  if("adjusted" %in% type & "regular" %in% type)
    return(list(regular = object$var,
         adjusted = object$var_adj))
  if(type == "adjusted")
    return(object$var_adj)
  if(type == "regular")
    return(object$var)
}

#' Residuals for frailty models
#'
#' @param object An \code{emfrail} object
#' @param type One of \code{cluster} or \code{individual}
#' @param ... Other arguments
#' @return A vector corresponding to the Martingale residuals, either for each cluster or for each individual (row of the data).
#'
#' @details For cluster \eqn{i}, individual \eqn{j} and observation row \eqn{k}, we write the cumulative hazard contribution as
#' \deqn{\Lambda_{ijk} = \exp(\beta^\top \mathbf{x}_{ijk}) \Lambda_{0, ijk}}
#' where \eqn{\Lambda_{0, ijk}} is the baseline cumulative hazard correspinding to the row \eqn{(i,j,k)}.
#'
#' When \code{type == "individual"}, the returned residuals are equal to \eqn{z_i \Lambda_{ijk}} where \eqn{z_i} is the estimated frailty in cluster \eqn{i}.
#' When \code{type == "cluster"}, the returned residuals are equal to \eqn{\sum_{j,k} \Lambda_{ijk}},
#' @export
#'
residuals.emfrail <- function(object, type = "group", ...) {
  object$residuals
}

#' @export
#' @method model.matrix emfrail
model.matrix.emfrail <- function(object, ...) {
  if(is.null(object$mm)) stop("emfrail must be called with model.matrix = TRUE in order to return the model matrix")
  object$mm
}


#' @export
#' @method model.frame emfrail
model.frame.emfrail <- function(formula, ...) {
  if(is.null(formula$mf)) stop("emfrail must be called with model = TRUE in order to return the model frame")
  formula$mf
}

#' Log-likelihood for \code{emfrail} fitted models
#'
#' @param object An \code{emfrail} object
#' @param ... Other arguments
#' @return An object of class \code{logLik} containing the marginal log-likelihood of the fitted model
#'
#' @details The formula for the likelihood can be found in the manual which accompanies the package. Note that a constant
#' is added. If we denote \eqn{\bar{n}} the total number of events and \eqn{\bar{n_i}} the total number of events at time point
#' \eqn{i}, for each time point where events are observed, then this is equal to
#' \deqn{\bar{n} - \sum_i \bar{n_i} \log \bar{n_i}.}
#' This is mostly because of compatibility, i.e. to match the log-likelihood given by the \code{survival} package.
#'
#' The \code{df} attribute of this object is equal to the number of regression coefficents plus 1.
#' In general, the number of degrees of freedom for a frailty model is an unclear concept. For the \code{coxph} frailty fits,
#' and in general for the shared frailty models fitted by penalized likelihood, the degrees of freedom is a number
#' that depends on the penalization. However, even in that case, there is no straight forward interpretation or use of this
#' quantity. The decision made here is because this would keep the likelihood ratio test for a covariate effect valid.
#' @export
logLik.emfrail <- function(object, ...) {
  res <- object$loglik[2]
  attr(res, "df") <- length(object$coefficients) + 1
  attr(res, "class") <- "logLik"
  res
}
