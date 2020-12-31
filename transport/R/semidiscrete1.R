#' Compute Semidiscrete Optimal Transport for Euclidean Distance Cost
#'
#' Computes the weight vector of the Apollonius diagram describing the semidiscrete
#' optimal transport plan for the Euclidean distance cost function and the associated
#' Wasserstein distance.
#'
#' @param source A matrix specifing the source measure.
#' @param target A three-column matrix specifing the target measure in the form
#'        x-coordinate, y-coordinate, mass.
#' @param xrange,yrange Vectors with two components defining the window on which 
#'        the source measure lives. Defaults to \eqn{[0,1] \times [0,1]}{[0,1] x [0,1]}.
#'        \code{source} is interpreted as an image of equally sized quadratic pixels
#'        on this window.
#' @param verbose Logical. Shall information about multiscale progress and L-BFGS return
#'        codes be printed?
#' @param reg A non-negative regularization parameter. It is usually not
#'        necessary to deviate from the default 0.
#'
#' @return A list describing the solution. The components are
#'   \item{weights}{A vector of length equal to the first dimension of \code{target}
#'         containing the weights for the Apollonius diagram discribing the
#'         optimal semidiscrete transport from source to target.}
#'   \item{wasserstein_dist}{The \eqn{L_1}{L1}-Wasserstein distance between source and target.}
#'   \item{ret_code}{A return code. Equal to 1 if everything is OK, since our code
#'         interrupts the usual lbfgs code. Other values can be converted to the
#'         corresponding return message by using \code{\link{ret_message}}.}
#' 
#' @note This function requires the Computational Geometry Algorithms Library (CGAL),
#'       available at \url{https://www.cgal.org}. Adapt the file src/Makevars according
#'       to the instructions given there and re-install from source.
#'       
#'       Internally the code from liblbfgs 1.10 by Naoaki Okazaki (2010) is used.
#'       See \url{http://www.chokkan.org/software/liblbfgs/}.
#'       
#'       A stand-alone version of the C++ code of this function is available
#'       at \url{https://github.com/valentin-hartmann-research/semi-discrete-transport}.  
#'
#' @author Valentin Hartmann \email{valentin.hartmann@epfl.ch} (stand-alone C++ code)\cr
#'         Dominic Schuhmacher \email{schuhmacher@math.uni-goettingen.de} (R-port)
#'
#' @references V. Hartmann and D. Schuhmacher (2017).
#'             Semi-discrete optimal transport --- the case p=1.
#'             Preprint \href{https://arxiv.org/abs/1706.07650}{arXiv:1706.07650}
#'             
#'             Menelaos Karavelas and Mariette Yvinec. 2D Apollonius Graphs 
#'             (Delaunay Graphs of Disks). In CGAL User and Reference Manual.
#'             CGAL Editorial Board, 4.12 edition, 2018
#'             
#'             Naoaki Okazaki (2010). libLBFGS: a library of Limited-memory
#'             Broyden-Fletcher-Goldfarb-Shanno (L-BFGS). Version 1.10
#'
#' @seealso \code{\link{ret_message}}, \code{\link{semidiscrete}}.
#'
#' @examples
#' \dontrun{
#' # the following function rotates a matrix m clockwise, so
#' # that image(rococlock(m)) has the same orientation as print(m):
#' roclock <- function(m) t(m)[, nrow(m):1]
#'
#' set.seed(30)
#' n <- 20
#' nu <- matrix(c(runif(2*n), rgamma(n,3,1)), n, 3)
#' pixelbdry <- seq(0,1,length=33)
#' image(pixelbdry, pixelbdry, roclock(random32a$mass), asp=1, col = grey(seq(0,1,length.out=32)))
#' points(nu[,1], nu[,2], pch=16, cex=sqrt(nu[,3])/2, col=2)
#'
#' res <- semidiscrete1(random32a$mass, nu)
#' plot_apollonius(nu[,1:2], res$weights, show_weights = FALSE, add = TRUE)
#' points(nu[,1], nu[,2], pch=16, cex=sqrt(nu[,3])/2, col=2)}
#'
#'
#' @export
semidiscrete1 <- function(source, target, xrange=c(0,1), yrange=c(0,1), verbose=FALSE, reg=0) {
  nocgal <- !as.logical(.Call('_transport_cgal_present', PACKAGE = 'transport'))
  if (nocgal) {
    stop("Computation of semidiscrete optimal transport for p=1 requires CGAL (the Computational
Geometry Algorithms Library). Install it from https://www.cgal.org/download.html and 
adapt the file src/Makevars of this package according to the instructions given there.
Then re-install 'transport' from source as usual.")
  }
  stopifnot(dim(target)[1] > 0 && dim(target)[2] == 3)
  xlen <- xrange[2] - xrange[1]
  ylen <- yrange[2] - yrange[1]
  maxlen <- max(xlen,ylen)
  rat_source <- dim(source)[2]/dim(source)[1]
  rat_range <- xlen/ylen
  if (!isTRUE(all.equal(rat_source,rat_range))) {
    stop("Source has aspect ratio ", rat_source, ", but ranges have aspect ratio ", rat_range)
  }
  stopifnot(all(xrange[1] <= target[,1], target[,1] <= xrange[2]))
  stopifnot(all(yrange[1] <= target[,2], target[,2] <= yrange[2]))

  if (dim(target)[1] == 1) {
    # there is no transport problem to solve and the C++ code would segfault
    psource <- source/sum(source)
    tx <- target[1,1]
    ty <- target[1,2]
    m <- dim(source)[2]
    n <- dim(source)[1]
    # rebuilding the refinement in opt_tranport.cpp / Source.cpp for consistency
    # of resulting Wasserstein distance
    refinement = ceiling(sqrt(1000 / m*n));
    mfin <- m * refinement
    nfin <- n * refinement
    sxgrid <- seq(xrange[1]+xlen/(2*mfin),xrange[2]-xlen/(2*mfin), length.out=mfin)
    sygrid <- seq(yrange[1]+ylen/(2*nfin),yrange[2]-ylen/(2*nfin), length.out=nfin)
    source_coord <- expand.grid(rev(sygrid),sxgrid)
    # the grouping matrix used also for the multiscale methods:
    grmat <- matrix(unlist(lapply(as.list(1:m), function(x) {rep(((x-1)*n+1):(x*n), 
                                        times = refinement, each = refinement)})), nfin, mfin)
    psource_refined <- matrix(psource[grmat]/refinement^2, refinement*n, refinement*m)
    cost <- sum(psource_refined * sqrt((source_coord[,2] - tx)^2 + (source_coord[,1] - ty)^2))
    res <- list(weights=0, wasserstein_dist=cost, ret_code=1)
    return(res)
  }
  
  target_norm <- target
  target_norm[,1] <- (target[,1]-xrange[1])/maxlen
  target_norm[,2] <- (target[,2]-yrange[1])/maxlen
  res <- semidiscrete_p1(source, target_norm, verbose, target_in_genpos = TRUE,
                             regularization_strength=reg, transportplan = matrix(0,1,1)) 
  res[[1]] <- maxlen*res[[1]] # weights are computed with distances scaled to [0,1]^2
                              # rescale to get back to distance scale
  res[[2]] <- maxlen*res[[2]] # same for resulting Wasserstein distance
  names(res) <- c("weights", "wasserstein_dist", "ret_code")
  return(res) 
}






#' Return Text Strings for lbfgs Return Codes
#' 
#' Given a vector of return codes, give back the corresponding vector of 
#' return strings from the lbfgs library. Nonexistant codes are ignored. 
#'
#' @param n The vector of return codes or \code{NULL} meaning the whole list
#'        shall be returned.
#'
#' @return A named character vector of the corresponding return strings.
#'  
#' @note Code 0 is ignored, since for technical reasons it is never returned by
#'       the function \code{\link{semidiscrete1}}.
#'
#' @author Dominic Schuhmacher \email{schuhmacher@math.uni-goettingen.de}
#'
#' @seealso \code{\link{semidiscrete1}}.
#'
#' @examples
#' ret_message()
#' ret_message(c(2,-1023,-1019))
#'
#' @keywords return code
#'
#' @export

ret_message <- function(n=NULL) {
  code <- list("1" = "LBFGS_STOP",
         "2" = "LBFGS_ALREADY_MINIMIZED",
         "-1024" = "LBFGSERR_UNKNOWNERROR",
         "-1023" = "LBFGSERR_LOGICERROR",
         "-1022" = "LBFGSERR_OUTOFMEMORY",
         "-1021" = "LBFGSERR_CANCELED",
         "-1020" = "LBFGSERR_INVALID_N",
         "-1019" = "LBFGSERR_INVALID_N_SSE",
         "-1018" = "LBFGSERR_INVALID_X_SSE",
         "-1017" = "LBFGSERR_INVALID_EPSILON",
         "-1016" = "LBFGSERR_INVALID_TESTPERIOD",
         "-1015" = "LBFGSERR_INVALID_DELTA",
         "-1014" = "LBFGSERR_INVALID_LINESEARCH",
         "-1013" = "LBFGSERR_INVALID_MINSTEP",
         "-1012" = "LBFGSERR_INVALID_MAXSTEP",
         "-1011" = "LBFGSERR_INVALID_FTOL",
         "-1010" = "LBFGSERR_INVALID_WOLFE",
         "-1009" = "LBFGSERR_INVALID_GTOL",
         "-1008" = "LBFGSERR_INVALID_XTOL",
         "-1007" = "LBFGSERR_INVALID_MAXLINESEARCH",
         "-1006" = "LBFGSERR_INVALID_ORTHANTWISE",
         "-1005" = "LBFGSERR_INVALID_ORTHANTWISE_START",
         "-1004" = "LBFGSERR_INVALID_ORTHANTWISE_END",
         "-1003" = "LBFGSERR_OUTOFINTERVAL",
         "-1002" = "LBFGSERR_INCORRECT_TMINMAX",
         "-1001" = "LBFGSERR_ROUNDING_ERROR",
         "-1000" = "LBFGSERR_MINIMUMSTEP",
         "-999" = "LBFGSERR_MAXIMUMSTEP",
         "-998" = "LBFGSERR_MAXIMUMLINESEARCH",
         "-997" = "LBFGSERR_MAXIMUMITERATION",
         "-996" = "LBFGSERR_WIDTHTOOSMALL",
         "-995" = "LBFGSERR_INVALIDPARAMETERS",
         "-994" = "LBFGSERR_INCREASEGRADIENT")
  if (is.null(n)) {
    subcode <- code
  } else {
    charn <- as.character(n)
    ind <- na.omit(match(charn,names(code)))
    if (length(ind) == 0) {
      warning("No existing return codes supplied")
      return(invisible(vector("character",0)))
    }
    subcode <- code[ind]
  }
  res <- as.character(subcode)
  names(res) <- names(subcode)
  return(res)
}


