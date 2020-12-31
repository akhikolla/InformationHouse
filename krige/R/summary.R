######################### SUMMARY FOR KRIGE OBJECT ############################
### LAST UPDATE: 11/08/2020; Le Bao

#' Summarize Fitted Kriging Model
#' 
#' Create a summary of a estimated model from \code{metropolis.krige}
#' 
#' @param object An \code{krige} object from the \code{metropolis.krige} function.
#' @param \dots Additional arguments passed to \code{summary} methods. Not supported for 
#'   \code{krige} object.
#'   
#' @details The function creates a summary of the model estimated by \code{metropolis.krige}.
#'   The output includes both the parameters and estimates of the model. 
#'   
#' @return A \code{summary.krige} list object.
#' 
#' @seealso \code{\link{as.mcmc.summary.krige}}
#' 
#' @examples
#' \dontrun{
#' # Summarize data
#' summary(ContrivedData)
#' 
#' # Initial OLS model
#' contrived.ols<-lm(y~x.1+x.2,data=ContrivedData)
#' # summary(contrived.ols)
#' 
#' # Set seed
#' set.seed(1241060320)
#' 
#' #For simple illustration, we set to few iterations.
#' #In this case, a 10,000-iteration run converges to the true parameters.
#' #If you have considerable time and hardware, delete the # on the next line.
#' #10,000 iterations took 39 min. with 8 GB RAM & a 1.5 GHz Quad-Core processor.
#' M <- 100
#' #M<-10000
#' 
#' contrived.run <- metropolis.krige(y ~ x.1 + x.2, coords = c("s.1","s.2"), 
#'    data = ContrivedData, n.iter = M, range.tol = 0.05)
#' 
#' # Summary
#' summary(contrived.run)
#' }
#' 
#' @importFrom stats quantile
#' @export
summary.krige <- function(object, ...) {
  if (!inherits(object, "krige")) stop("The input object is not a 'krige' object.")
  cl <- match.call()
  formula <- object$formula
  coords <- object$coords
  n.iter <- object$n.iter
  n.burnin <- object$n.burnin
  mcmc.mat <- object$mcmc.mat
  out <- apply(object$mcmc.mat, 2, quantile, c(0.5,0.0005,0.005,0.025,0.05,0.1,0.9,0.95,0.975,0.995,0.9995))
  estimates <- t(out); colnames(estimates)[1] <- "median"
  acceptance.rate <- object$acceptance.rate
  sumKrige <- list(call = cl, formula = formula, coords = coords, n.iter = n.iter, 
                   n.burnin = n.burnin, mcmc.mat = mcmc.mat, estimates = estimates, 
                   acceptance.rate = acceptance.rate)
  class(sumKrige) <- "summary.krige"
  sumKrige
} 


################ PRINT: PRINT RESUTLS OF metropolis.krige() ####################
### LAST UPDATED: 11/08/2020; Le Bao

#' @export
print.krige <- function(x, ..., accept.rate.warning, digits = max(3L, getOption("digits") - 3L)) {
  if (!inherits(x, "krige")) stop("The input object is not a 'krige' object.")
  
  if (missing(accept.rate.warning)) {
    if ("accept.rate.warning" %in% names(x$call)) {
      accept.rate.warning <- x$call$accept.rate.warning} else {
        accept.rate.warning <- TRUE
      } }
  
  cat("Call: ", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")
  cat("Acceptance rates:", "beta:", paste0(round(x$acceptance.rate[["beta.rate"]],digits)*100, "%"),
      "tau2:", paste0(round(x$acceptance.rate[["tau2.rate"]],digits)*100, "%"),
      "phi:", paste0(round(x$acceptance.rate[["phi.rate"]],digits)*100, "%"),
      "sigma:", paste0(round(x$acceptance.rate[["sigma2.rate"]],digits)*100, "%."))
  
  if (accept.rate.warning==TRUE) ar.check(x$acceptance.rate)
  invisible(x)
}

################  PRINT: PRINT RESUTLS OF summary.krige() ######################
### LAST UPDATED: 11/08/2020; Le Bao

#' @export
print.summary.krige <- function(x, digits=max(3, getOption("digits") - 3),..., accept.rate.warning) {
  if (!inherits(x, "summary.krige")) stop("The input object is not a 'summary.krige' object.")
  if (missing(accept.rate.warning)) {
    if ("accept.rate.warning" %in% names(x$call)) {
      accept.rate.warning <- x$call$accept.rate.warning} else {
        accept.rate.warning <- TRUE
      } }
  cat("Formula:", paste(deparse(x$formula), sep = "\n", collapse = "\n"), 
      "\n")
  cat("\nCoordinates: ", paste0('"', x$coords[1],'"'), ",", paste0('"', x$coords[2],'"'), 
      "\n", sep = "", collapse = "")
  cat("\nNumber of iterations: ", x$n.iter, " (burnin: ", x$n.burnin, ")", "\n", sep = "")
  cat("\nEstimates:\n")
  sumTab <- round(apply(x$mcmc.mat, 2, quantile, c(0.5,0.025,0.05,0.95,0.975)), digits)
  outTab <- t(sumTab); colnames(outTab)[1] <- "median"
  print(outTab)
  cat("---\n")
  cat("Acceptance rates:", "beta:", paste0(round(x$acceptance.rate[["beta.rate"]],3)*100, "%"),
      "tau2:", paste0(round(x$acceptance.rate[["tau2.rate"]],3)*100, "%"),
      "phi:", paste0(round(x$acceptance.rate[["phi.rate"]],3)*100, "%"),
      "sigma:", paste0(round(x$acceptance.rate[["sigma2.rate"]],3)*100, "%."))
  if (accept.rate.warning==TRUE) ar.check(x$acceptance.rate)
}

############################## BURN IN FUNCTION ################################

#' Discard Burn-in Period of Kriging Model
#' 
#' Discard burn-in period of a estimated model from \code{metropolis.krige}
#' 
#' @param object An \code{krige} object from the \code{metropolis.krige} function.
#' @param n.burnin The number of burnin iterations. Defaults to half of the iterations.
#'   
#' @details The function discard the burn-in period from the results of \code{metropolis.krige}.
#'   It is generally used for discarding burn-in for \code{krige} object that keeps 
#'   all the iterations.
#' 
#' @importFrom utils tail
#' @export
burnin <- function(object, n.burnin) UseMethod("burnin")

#' @rdname burnin
#' @method burnin krige
#' @export
burnin.krige <- function(object, n.burnin=object$n.iter/2){
  if (!inherits(object, "krige")) stop("The input object is not a 'krige' object.")
  if (object$n.burnin != 0) warning("'n.burnin' is already specified in 'krige' object. Additional burn-in iterations will be discarded.")
  if (n.burnin < 0) {
    n.burnin <- 0; warning("The burn-in period is negative. 'n.burnin = 0' is used.")
  } else if (n.burnin >= object$n.iter) {stop("The number of iterations is less than the burn-in period.")}
  if (n.burnin%%1 != 0) {
    n.burnin <- round(n.burnin) 
    warning("The number of burn-in is not an integer. 'round(n.burnin)' is used.")
  }
  object$call$n.burnin <- object$n.burnin + n.burnin
  object$mcmc.mat <- tail(object$mcmc.mat, nrow(object$mcmc.mat) - n.burnin)
  object$n.burnin <- object$n.burnin + n.burnin
  object
}

#' @rdname burnin
#' @method burnin matrix
#' @export
burnin.matrix <- function(object, n.burnin=nrow(object)/2){
  if (!inherits(object, "matrix")) stop("The input object is not a 'matrix' object.")
  if (n.burnin < 0) {
    n.burnin <- 0; warning("The burn-in period is negative. 'n.burnin = 0' is used.")
  } else if (n.burnin >= nrow(object)) {stop("The number of iterations is less than the burn-in period.")}
  if (n.burnin%%1 != 0) {
    n.burnin <- round(n.burnin) 
    warning("The number of burn-in is not an integer. 'round(n.burnin)' is used.")
  }
  object <- tail(object, nrow(object) - n.burnin)
  object
}