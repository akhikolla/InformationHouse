#' Compute Optimal Transport (Cost/Plan) Using the Multiscale Shielding Method
#' 
#' Runs the multiscale version of the Shielding Method (a.k.a. Short Cut Method) for computing the optimal transport
#' (cost/plan) on a rectangular grid in \eqn{d} dimensions for the squared Euclidean distance as cost function.
#'
#' If \code{a} and \code{b} do not have the same sum, they are normalized to sum 1 \emph{before}
#' \code{flood} and \code{measureScale} transformations are applied.
#'
#' @section Use of CPLEX:
#' For larger problems (thousands of grid points) there are considerable speed improvements when \code{shielding}
#' can use the CPLEX numerical solver for the underlying constrained optimization problems.
#' If a local installation of CPLEX is available, the transport package can be linked against it during installation.
#' See the file src/Makevars in the source package for instructions.
#'
#' @param a,b arrays with \eqn{d} coordinates representing source and target measure, respectively.
#'        The entries must be all positive.  
#' @param nscales the number of scales generated in the multiscale algorithm.
#' @param startscale the first scale on which the problem is solved.
#' @param flood a real number. If positive, take the maximum of entry and \code{flood} for each
#'        entry of \code{a} and \code{b}.
#' @param measureScale the required precision for the entries. Computations are performed on
#'        \code{round(a/measureScale)} and the same for \code{b} using integer arithmetics. 
#' @param verbose logical. Toggles output to the console about the progress of the algorithm.        
#' @param basisKeep,basisRefine internal use only.                               
#' 
#' @return A list of components
#' \item{err}{error code. 0 if everything is ok.} 
#' \item{a_used,b_used}{the vectorized arrays that were actually used by the algorithm. \code{a}, \code{b} after
#'                 applying \code{flood} and \code{measureScale}.}  
#' \item{coupling}{a vectorized coupling describing the optimal transport from a_used to b_used}
#' \item{basis}{a matrix with two columns describing the basis obtained for the optimal transport}
#' \item{u,v}{vectors of optimal values in the dual problem}
#'
#' @seealso \code{\link{transport}}, which calls this function if appropriate.
#'
#' @references B. Schmitzer (2016). A sparse multiscale algorithm for dense optimal transport. J. Math. Imaging Vision 56(2), 238--259. \href{https://arxiv.org/abs/1510.05466}{https://arxiv.org/abs/1510.05466}
#'
#' @author Bernhard Schmitzer \email{schmitzer@uni-muenster.de} and\cr
#'         Dominic Schuhmacher \email{dschuhm1@uni-goettingen.de}\cr
#'         (based on C++ code by Bernhard Schmitzer)
#'
#' @examples
#' \dontrun{
#' shielding(random64a$mass,random64b$mass,nscales=6,measureScale=1) }
#'
#' @export
shielding <- function(a,b,nscales=2,startscale=1,flood=0,measureScale=1e-6,verbose=FALSE,basisKeep=1,basisRefine=1) {
  # measureScale is precision for the measures as they are entered, 1e-6 up to sixth digit behind decimal point
  # the user has to take care that the sum of all masses in a and b separately if divided by measure scale
  # is smaller equal the largest representable integer .Machine$integer.max
  # This should be changed: it can be done in such a way that only the max mass divided by measures scale has to be
  # smaller equal the sum
  stopifnot(is.array(a) && is.array(b))
  stopifnot(all(dim(a) == dim(b)))
  if (sum(a) != sum(b)) {
    a <- a/sum(a)
    b <- b/sum(b)	
  }
  stopifnot(startscale < nscales)

  nocplex <- !as.logical(.Call('_transport_cplex_present', PACKAGE = 'transport'))
  if (nocplex && !isTRUE(getOption("transport-CPLEX_no_warn"))) {
     message("If you have CPLEX available, the computation of optimal transport for the current data (with p=2) can be accelerated considerably by linking against it.
To do so, adapt the file src/Makevars in the package source of 'transport' according to the instructions given there. Then install from source as usual.
This message will appear only once per session. To turn it off completely say 'options(\"transport-CPLEX_no_warn\" = TRUE)' on startup (e.g. in a .Rprofile file).")
    options("transport-CPLEX_no_warn" = TRUE)
  }
  
  if (flood > 0) {	  	
    a <- pmax(a,flood)
    b <- pmax(b,flood)
  }

  m <- length(a)
  n <- length(b)
  assignment <- matrix(0,m,n)
  u <- rep(0,m)
  v <- rep(0,n)

  res <- .Call('_transport_SolveHierarchicalTransport', PACKAGE = 'transport',
               x=as.numeric(a), y=as.numeric(b), xydepth=length(dim(a)), xydimensions=dim(a),
               compdepth=as.integer(nscales), measureScaleVecPre=as.double(measureScale),
               keepBasisVecPre=as.integer(basisKeep), refineBasisVecPre=as.integer(basisRefine),
               layerCoarsestVecPre=as.integer(startscale), verboseVecPre=verbose,
               assignment=as.integer(assignment), u=as.double(u), v=as.double(v))
  if (res[[1]] != 0) 
    warning("SolveHierarchicalTransport (method='shielding') returned an error. The result may be incorrect.")
  names(res) <- c("err","a_used","b_used","coupling","basis","u","v")
  return(res)
}


costp2 <- function(m,n) {
  gg <- expand.grid(1:m,1:n)
  dd <- as.matrix(dist(gg))^2
  return(dd)
}


res2matrix <- function(res,m,n) {
  k <- dim(res)[1]
  resmat <- matrix(0,m*n,m*n)
  for (i in 1:k) {
  	resmat[res[i,1],res[i,2]] <- res[i,3]
  }
  return(resmat)
}
