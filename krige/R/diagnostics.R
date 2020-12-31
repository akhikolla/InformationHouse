#################### S3 METHOD FOR CREATING MCMC OBJECT #######################
### LAST UPDATE: 11/08/2020; Le Bao

#' Convert \code{krige} object to an \code{mcmc} object
#' 
#' Convert MCMC matrix of posterior samples for use with the \pkg{coda} package
#' 
#' @param x An \code{krige} or \code{summary.krige} object.
#' @param start The iteration number of the first observation.
#' @param end The iteration number of the last observation.
#' @param thin The thinning interval between consecutive observations.
#' @param \dots Additional arguments to be passed to \code{mcmc()} methods of \code{coda}
#'   package.
#'   
#' @details The function converts a \code{krige} output object to a Markov Chain 
#'   Monte Carlo (mcmc) object used in \code{coda} as well as a variety of MCMC 
#'   packages. It extracts the MCMC matrix of posterior samples from the output 
#'   of \code{metropolis.krige} for further use with other MCMC packages and functions. 
#'   
#' @return A \code{mcmc} object.
#' 
#' @seealso \code{\link[coda:as.mcmc]{coda::as.mcmc()}}
#' 
#' @examples
#' \dontrun{
#' # Summarize Data
#' summary(ContrivedData)
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
#'                                   data = ContrivedData, n.iter = M, n.burnin = 20,
#'                                   range.tol = 0.05)
#'                                   
#' # Convert to mcmc object
#' mcmc.contrived.run <- as.mcmc(contrived.run)
#' #mcmc.contrived.run <- as.mcmc(summary(contrived.run))
#' 
#' # Diagnostics using MCMC packages
#' coda::raftery.diag(mcmc.contrived.run)
#' # superdiag::superdiag(mcmc.contrived.run) #NOT WORKING YET
#' }
#' 
#' @method as.mcmc krige
#' @importFrom coda mcmc as.mcmc
#' @export
as.mcmc.krige <- function(x, start = 1, end = x$n.iter, thin = 1, ...){
  if (!inherits(x, "krige")) stop("The input object is not a 'krige' object.")
  if (! requireNamespace("coda", quietly = TRUE)) {
    stop("There is no package called 'coda'. Install or load 'coda' package first.")
  }
  mcmc(data = x$mcmc.mat, start = start, end = end, thin = thin, ...)
}

#' @rdname as.mcmc.krige
#' @method as.mcmc summary.krige
#' @export
as.mcmc.summary.krige <- function(x, start=1, end = x$n.iter, thin = 1, ...){
  if (!inherits(x, "summary.krige")) stop("The input object is not a 'summary.krige' object.")
  if (! requireNamespace("coda", quietly = TRUE)) {
    stop("There is no package called 'coda'. Install or load 'coda' package first.")
  }
  mcmc(data = x$mcmc.mat, start = start, end = end, thin = thin, ...)
} 

######################### GEWEKE DIAGNOSTIC FUNCTION ###########################
# Note: A simplified version
# LAST UPDATE: 11/15/2020

#' Geweke Diagnostic for MCMC
#' 
#' Conducts a Geweke convergence diagnostic on MCMC iterations.
#' 
#' @param object An matrix or \code{krige}/\code{summary.krige} object for which 
#' a Geweke diagnostic is desired
#' @param early.prop Proportion of iterations to use from the start of the chain.
#' @param late.prop Proportion of iterations to use from the end of the chain.
#' @param precision Number of digits of test statistics and p-values to print.
#' 
#' @details This is a generic function currently works with \code{matrix}, \code{krige}, 
#'   and \code{summary.krige} objects. It is a simplified version of the Geweke 
#'   test for use with this package.
#'   
#'   Geweke's (1992) test for nonconvergence of a MCMC chain is to conduct a 
#'   difference-of-means test that compares the mean early in the chain to the mean 
#'   late in the chain. If the means are significantly different from each other, 
#'   then this is evidence that the chain has not converged. The difference-of-means 
#'   test is a simple z-ratio, though the standard error is estimated using the 
#'   spectral density at zero to account for autocorrelation in the chain.
#'   
#' @return A \code{matrix} in which the first row consists of z-scores for tests 
#'   of equal means for the first and last parts of the chain. The second row 
#'   consists of the corresponding p-values. Each column of the matrix represents 
#'   another parameter. For each column, a significant result is evidence that the 
#'   chain has not converged for that parameter. Thus, a non-significant result 
#'   is desired.
#' 
#' @references 
#'   John Geweke. 1992. "Evaluating the Accuracy of Sampling-Based Approaches to 
#'   the Calculation of Posterior Moments." In \emph{Bayesian Statistics 4}, ed. 
#'   J.M. Bernardo, J.O. Berger, A.P. Dawid, and A.F.M. Smith. Oxford: Clarendon 
#'   Press.
#' 
#' @seealso \code{\link{geweke.krige}}, \code{\link{geweke.summary.krige}}
#' 
#' @rdname geweke
#' 
#' @examples
#' \dontrun{
#' # Load Data
#' data(ContrivedData)
#' 
#' # Set seed
#' set.seed(1241060320)
#' 
#' M <- 100
#' 
#' contrived.run <- metropolis.krige(y ~ x.1 + x.2, coords = c("s.1","s.2"), 
#'   data = ContrivedData, n.iter = M, range.tol = 0.05)
#' 
#' geweke(contrived.run, early.prop=0.5)
#' geweke(summary(contrived.run), early.prop=0.5)
#' geweke(contrived.run$mcmc.mat, early.prop=0.5)
#' # Note that the default proportions in the geweke() is more typical for longer run.
#' }
#' 
#' @importFrom stats ar pnorm
#' @export
geweke<-function(object,early.prop=.1,late.prop=.5,precision=4){
  UseMethod("geweke")
}

#' @rdname geweke
#' @export
geweke.krige<-function(object,early.prop=.1,late.prop=.5,precision=4){ 
  if (!inherits(object, "krige")) stop("The input object is not a 'krige' object.")
  X <- as.matrix(object$mcmc.mat)
  geweke.default(X, early.prop=.1,late.prop=.5,precision=4)
}

#' @rdname geweke
#' @export
geweke.summary.krige<-function(object,early.prop=.1,late.prop=.5,precision=4){
  if (!inherits(object, "summary.krige")) stop("The input object is not a 'summary.krige' object.")
  X <- as.matrix(object$mcmc.mat)
  geweke.default(X,early.prop=.1,late.prop=.5,precision=4)
}

#' @rdname geweke
#' @export
geweke.default<-function(object,early.prop=.1,late.prop=.5,precision=4){ 
  if (!is.matrix(object)) {
    object <- as.matrix(object); warning("Input object was not matrix class. Coerced to a matrix.")}  
  if(early.prop<0 || early.prop>1) stop("early.prop must be between 0 and 1.")
  if(late.prop<0 || late.prop>1) stop("late.prop must be between 0 and 1.")
  if(early.prop+late.prop>1) stop("The sum of early.prop+late.prop must be less than 1.")
  if(precision<=2) stop("precision must be 2 or larger.")
  X <- object
  iter<-dim(X)[1]
  K<-dim(X)[2]
  early<-c(1:ceiling(1+early.prop*(iter-1)))
  late<-c(floor(iter-late.prop*(iter-1)):iter)
  n.early<-length(early)
  n.late<-length(late)
  geweke.vec<-rep(NA,K)
  for(k in 1:K){
    early.ar<-ar(X[early,k],aic=T)
    early.var<-early.ar$var.pred/(1-sum(early.ar$ar))^2
    late.ar<-ar(X[late,k],aic=T)
    late.var<-late.ar$var.pred/(1-sum(late.ar$ar))^2
    geweke.vec[k]<-(mean(X[early,k])-mean(X[late,k]))/sqrt((early.var/n.early)+(late.var/n.late))
  }
  geweke.mat<-rbind(round(geweke.vec,precision),round(2*(1-pnorm(abs(geweke.vec))),precision))
  colnames(geweke.mat)<-colnames(X)
  row.names(geweke.mat)<-c("z-ratio","p-value")
  return(geweke.mat)
}

################# HEIDELBERGER AND WELCH DIAGNOSTIC FUNCTION ###################
# Note: A simplified version adapted from the heidel.diag() function in the coda package
# LAST UPDATE: 11/15/2020

#' Heidelberger and Welch Diagnostic for MCMC
#' 
#' Conducts a Heidelberger and Welch convergence diagnostic on MCMC iterations.
#' 
#' @param object An matrix or \code{krige}/\code{summary.krige} object for which 
#' a Heidelberger and Welch diagnostic is desired
#' @param pvalue Alpha level for significance tests. Defaults to 0.05.
#'   
#' @details This is a generic function currently works with \code{matrix}, \code{krige}, 
#'   and \code{summary.krige} objects. It is a simplified version of the Heidelberger 
#'   and Welch test for use with this package.
#'   
#'    This is an adaptation of a function in Plummer et al.'s \code{coda} package. 
#'    Heidelberger and Welch's (1993) test for nonconvergence. This version of the 
#'    diagnostic only reports a Cramer-von Mises test and its corresponding p-value 
#'    to determine if the chain is weakly stationary with comparisons of early 
#'    portions of the chain to the end of the chain.
#'    
#' @return A \code{matrix} in which the first row consists of the values of the 
#'   Cramer-von Mises test statistic for each parameter, and the second row consists 
#'   of the corresponding p-values. Each column of the matrix represents another 
#'   parameter of interest. A significant result serves as evidence of nonconvergence, 
#'   so non-significant results are desired.
#' 
#' @references 
#' Philip Heidelberger and Peter D. Welch. 1993. "Simulation Run Length Control 
#' in the Presence of an Initial Transient." \emph{Operations Research} 31:1109-1144.
#' 
#' Martyn Plummer, Nicky Best, Kate Cowles and Karen Vines. 2006. "CODA: Convergence 
#' Diagnosis and Output Analysis for MCMC." \emph{R News} 6:7-11.
#' 
#' @seealso \code{\link{heidel.welch.krige}}, \code{\link{heidel.welch.summary.krige}}, 
#' \code{\link{geweke}}
#' 
#' @examples
#' \dontrun{
#' # Load Data
#' data(ContrivedData)
#' 
#' # Set seed
#' set.seed(1241060320)
#' 
#' M <- 100
#' 
#' contrived.run <- metropolis.krige(y ~ x.1 + x.2, coords = c("s.1","s.2"), 
#'   data = ContrivedData, n.iter = M, n.burnin = 20, range.tol = 0.05)
#' 
#' heidel.welch(contrived.run)
#' heidel.welch(summary(contrived.run))
#' heidel.welch(contrived.run$mcmc.mat)
#' }
#' 
#' @importFrom stats start end window
#' @export

heidel.welch<-function(object, pvalue){
  UseMethod("heidel.welch")
}

#' @rdname heidel.welch
#' @export

heidel.welch.krige<-function(object, pvalue=0.05){
  if (!inherits(object, "krige")) stop("The input object is not a 'krige' object.")
  heidel.welch.default(as.matrix(object$mcmc.mat), pvalue = pvalue)
}

#' @rdname heidel.welch
#' @export
heidel.welch.summary.krige<-function(object, pvalue=0.05){
  if (!inherits(object, "summary.krige")) stop("The input object is not a 'summary.krige' object.")
  heidel.welch.default(as.matrix(object$mcmc.mat), pvalue = pvalue)
}  

#' @rdname heidel.welch
#' @export

heidel.welch.default<-function(object, pvalue=0.05){
  if (!is.matrix(object)) {
    object <- as.matrix(object); warning("Input object was not matrix class. Coerced to a matrix.")}
  X <- object
  HW.mat <- matrix(0,ncol=ncol(X),nrow=2)
  dimnames(HW.mat) <- list(c("CramerVonMises", "p-value"), colnames(X))
  pcramer<-function(q, eps = 1e-05) {
    log.eps <- log(eps)
    y <- matrix(0, nrow = 4, ncol = length(q))
    for (k in 0:3) {
      z <- gamma(k + 0.5) * sqrt(4 * k + 1)/(gamma(k + 1)*pi^(3/2) * sqrt(q))
      u <- (4 * k + 1)^2/(16 * q)
      y[k + 1, ] <- ifelse(u > -log.eps, 0, z * exp(-u) * besselK(x = u,nu = 1/4))
    }
    return(apply(y, 2, sum))
  }
  
  for (j in 1:ncol(X)) {
    start.vec <- seq(from = start(X)[1], to = end(X)[1]/2, by = nrow(X)/10)
    Y <- X[, j, drop = TRUE]
    n1 <- length(Y)
    first.ar<-ar(window(Y, start = end(Y)/2),aic=T)
    S0<-first.ar$var.pred/(1-sum(first.ar$ar))^2
    converged <- FALSE
    for (i in seq(along = start.vec)) {
      Y <- window(Y, start = start.vec[i])
      n <- length(Y)
      ybar <- mean(Y)
      B <- cumsum(Y)-ybar*(1:n)
      Bsq <- (B*B)/(n*S0)
      I <- sum(Bsq)/n
      if (converged <- !is.na(I) && pcramer(I) < 1 - pvalue) 
        break
    }
    HW.mat[,j] <- c(I, 1 - pcramer(I))
  }
  return(HW.mat)
}


