################################# MCMC SAMPLES #################################
### LAST UPDATED: 11/08/2020; Le Bao

#' Extract MCMC Samples
#' 
#' Extract MCMC samples estimated by \code{metropolis.krige()}
#' 
#' @param object A \code{krige}or \code{summary.krige} object from the \code{metropolis.krige} function.
#' @param as.matrix Logical values indicating if the output format should be a matrix. Defaults to \code{TRUE}.
#' @param as.data.frame Logical values indicating if the output format should be a 
#'   data.frame. Defaults to \code{FALSE}.
#' @param x A \code{krige} or \code{summary.krige} object for \code{as.matrix} and \code{as.data.frame} methods.
#' @param \dots Additional arguments passed to \code{as.matrix} or \code{as.data.frame} methods. 
#'   
#' @details The function extracts the MCMC samples from the a \code{krige}or \code{summary.krige} 
#'   object from the \code{metropolis.krige} function. Users can define the output by using \code{as.matrix}
#'   or \code{as.data.frame}.
#'   
#' @return A \code{summary.krige} list object.
#' 
#' @seealso \code{\link{as.mcmc.krige}}
#' 
#' @examples
#' \dontrun{
#' # Summarize Data
#' summary(ContrivedData)
#' 
#' # Initial OLS model
#' contrived.ols<-lm(y~x.1+x.2,data=ContrivedData)
#' # summary(contrived.ols)
#' 
#' # Set seed
#' set.seed(1241060320)
#' 
#' M <- 100
#' #M<-10000
#' 
#' contrived.run <- metropolis.krige(y ~ x.1 + x.2, coords = c("s.1","s.2"), 
#'    data = ContrivedData, n.iter = M, n.burnin = 20, range.tol = 0.05)
#'    
#' contrived.run.mat <- mcmc.samples(contrived.run)
#' 
#' ### Alternatively, use generic methods
#' contrived.run.mat <- as.matrix(contrived.run)
#' contrived.run.df <- as.data.frame(contrived.run)
#' }
#' 
#' @importFrom utils tail
#' @export

mcmc.samples <- function(object, as.matrix, as.data.frame, ...) {
  UseMethod("mcmc.samples")
}

#' @rdname mcmc.samples
#' @export
mcmc.samples.krige <- function(object, as.matrix=!as.data.frame, as.data.frame=FALSE, ...) {
  if (as.matrix && as.data.frame) {
    stop("Cannot use 'as.matrix' and 'as.data.frame' at the same time.")
  }
  if (!inherits(object, "krige")) stop("The input object is not a 'krige' object.")
  out <- object$mcmc.mat
  if (as.data.frame) {
    out <- as.data.frame(out, ...)
  } else {
    out <- as.matrix(out, ...)
  }
  out
}

#' @rdname mcmc.samples
#' @export
mcmc.samples.summary.krige <- function(object, as.matrix=!as.data.frame, as.data.frame=FALSE, ...) {
  if (as.matrix && as.data.frame) {
    stop("Cannot use 'as.matrix' and 'as.data.frame' at the same time.")
  }
  if (!inherits(object, "summary.krige")) stop("The input object is not a 'summary.krige' object.")
  out <- object$mcmc.mat
  if (as.data.frame) {
    out <- as.data.frame(out, ...)
  } else {
    out <- as.matrix(out, ...)
  }
  out
}

#' @rdname mcmc.samples
#' @export
as.matrix.krige <- function(x, ...) {
  mcmc.samples(x, as.matrix=TRUE, ...)
}


#' @rdname mcmc.samples
#' @export
as.matrix.summary.krige <- function(x, ...) {
  mcmc.samples(x, as.matrix=TRUE, ...)
}

#' @rdname mcmc.samples
#' @export
as.data.frame.krige <- function(x, ...) {
  mcmc.samples(x, as.data.frame=TRUE, ...)
}

#' @rdname mcmc.samples
#' @export
as.data.frame.summary.krige <- function(x, ...) {
  mcmc.samples(x, as.data.frame=TRUE, ...)
}