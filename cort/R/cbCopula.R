#' @include generics.R empiricalCopula.R
NULL

.cbCopula = setClass(Class = "cbCopula", contains = "empiricalCopula",
  slots = c(m = "numeric"), validity = function(object) {
   errors <- c()
   if (any(sapply(object@m,function(m){(nrow(object@data)%%m) != 0}))) {
     errors <- c(errors, "m should divide the number of row")
   }
   if(length(object@m) != ncol(object@data)){
     errors <- c(errors, "lengths of m parameter shoudl be the same as the number of columns of data")
   }
   if (length(errors) == 0)
     TRUE else errors
  })

#' cbCopula contructor
#'
#' The cbCopula class computes a checkerboard copula with a given checkerboard parameter \eqn{m}, as described by A. Cuberos, E. Masiello and V. Maume-Deschamps (2019).
#' Assymptotics for this model are given by C. Genest, J. Neslehova and R. bruno (2017). The construction of this copula model is as follows :
#'
#' Start from a dataset with \eqn{n} i.i.d observation of a \eqn{d}-dimensional copula (or pseudo-observations), and a checkerboard parameter \eqn{m},dividing \eqn{n}.
#'
#' Consider the ensemble of multi-indexes \eqn{I = \{i = (i_1,..,i_d) \subset \{1,...,m \}^d\}} which indexes the boxes :
#'
#' \deqn{B_{i} = \left]\frac{i-1}{m},\frac{i}{m}\right]}
#'
# 'partitioning the space \eqn{I^d = (0,1)^d}.
#'
#' Let now \eqn{\lambda} be the dimension-unspecific lebesgue measure on any power of \eqn{R}, that is :
#'
#' \deqn{\forall d \in N, \forall x,y \in R^p, \lambda(\left(x,y\right)) = \prod\limits_{p=1}^{d} (y_i - x_i)}
#'
#' Let furthermore \eqn{\mu} and \eqn{\hat{\mu}} be respectively the true copula measure of the sample at hand and the classical Deheuvels empirical copula, that is :
#'
#' - For \eqn{n} i.i.d observation of the copula of dimension \eqn{d}, let \eqn{\forall i \in \{1,...,d\}, \, R_i^1,...,R_i^d} be the marginal ranks for the variable \eqn{i}.
#' - \eqn{\forall x \in I^d} let \eqn{\hat{\mu}((0,x)) = \frac{1}{n} \sum\limits_{k=1}^n I_{R_1^k\le x_1,...,R_d^k\le x_d}}
#'
#'
#' The checkerboard copula, \eqn{C}, and the empirical checkerboard copula, \eqn{\hat{C}}, are then defined by the following :
#'
#' \deqn{\forall x \in (0,1)^d, C(x) = \sum\limits_{i\in I} {m^d \mu(B_{i}) \lambda((0,x)\cap B_{i})}}
#'
#' Where \eqn{m^d = \lambda(B_{i})}.
#'
#' This copula is a special form of patchwork copulas, see F. Durante, J. Fernández Sánchez and C. Sempi (2013) and F. Durante, J. Fernández Sánchez, J. Quesada-Molina and M. Ubeda-Flores (2015).
#' The estimator has the good property of always being a copula.
#'
#' @param x the data to be used
#' @param m checkerboard parameters
#' @param pseudo Boolean, defaults to `FALSE`. Set to `TRUE` if you are already
#'  providing pseudo data into the `x` argument.
#'
#' @details The checkerboard copula is a kind of patchwork copula that only uses independent copula as fill-in, only where there are values on the empirical data provided.
#' To create such a copula, you should provide data and checkerboard parameters (depending on the dimension of the data).
#'
#' @name cbCopula-Class
#' @title Checkerboard copulas
#' @rdname cbCopula-Class
#'
#' @return An instance of the `cbCopula` S4 class. The object represent the fitted copula and can be used through several methods to query classical (r/d/p/v)Copula methods, etc.
#' @export
#'
#' @references
#' \insertRef{cuberos2019}{cort}
#'
#' \insertRef{genest2017}{cort}
#'
#' \insertRef{durante2013}{cort}
#'
#' \insertRef{durante2015}{cort}
#'
cbCopula = function(x, m = rep(nrow(x),ncol(x)), pseudo = FALSE) {
  if (missing(x)) { stop("The argument x must be provided") }
  if(ncol(x) == 0){ stop("you are providing a matrix equal to NULL") }
  if(nrow(x) == 0){ stop("You provided no data !") }
  if (!pseudo) {
    x <- apply(x, 2, rank, na.last = "keep")/(nrow(x) + 1)
  }
  if(length(m) == 1){m = rep(m,ncol(x))}
  if(length(m) != ncol(x)){stop("You should provide m values same lengths as the number of columns in data.")}
  return(.cbCopula(data = as.matrix(x), dim=ncol(x),m = m))
}
setMethod(f = "show",    signature = c(object = "cbCopula"),                definition = function(object)    {
  cat("This is a cbCopula , with : \n", "  dim =", dim(object), "\n   n =",
      nrow(object@data), "\n   m =", object@m, "\n")
  cat("The variables names are : ", colnames(object@data))
})


#' @describeIn rCopula-methods Method for the cbCopula
setMethod(f = "rCopula", signature = c(n = "numeric", copula = "cbCopula"), definition = function(n, copula) {

  # if n=0, return a 0xdim matrix :
  if (n == 0) { return(matrix(0, nrow = 0, ncol = copula@dim)) }
  x = copula@data
  d = copula@dim
  m = copula@m
  rows <- resample(x = 1:nrow(x), size = n, replace = TRUE)
  seuil_inf <- boxes_from_points(x,m)[rows,,drop=FALSE]
  size = matrix(rep(1/m,each = n),nrow=n,ncol=d)
  rng <- matrix(runif(d * n), nrow = n, ncol = d)
  return(seuil_inf + rng * size)
})

#' @describeIn pCopula-methods Method for the cbCopula
setMethod(f = "pCopula", signature = c(u = "matrix",  copula = "cbCopula"), definition = function(u, copula) {
  # The implementation is in Rcpp
  pcbCopula(u,copula@data,copula@m)
})

#' @describeIn dCopula-methods Method for the cbCopula
setMethod(f = "dCopula", signature = c(u = "matrix",  copula="cbCopula"),   definition = function(u, copula) {
  # The implementaiton is in Rcpp
  dcbCopula(u,copula@data,copula@m)
})




























