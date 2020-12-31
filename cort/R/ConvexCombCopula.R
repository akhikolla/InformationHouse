#' @include generics.R
NULL


############################### ConvexComCopula class ####### it is really a mixture of copulas.
.ConvexCombCopula = setClass(Class = "ConvexCombCopula",
  slots = c(copulas = "list", alpha = "numeric",dim="numeric"), validity = function(object) {
    errors <- c()
    if (length(object@copulas) != length(object@alpha)) {
      errors <- c(errors, "the weights parameter alpha must have same length as the copulas list")
    }
    if (!all(object@alpha >= 0)) {
      errors <- c(errors, "weights should be positives")
    }
    if(abs(sum(object@alpha)-1)>10^(-7)){
      errors <- c(errors,"weights should add up to 1.")
    }
    if (!(length(unique(sapply(object@copulas, dim))) == 1)) {
      errors <- c(errors, "all copulas must have same dimension")
    }
    if (length(errors) == 0)
      TRUE else errors
  })

#' ConvexCombCopula class
#'
#' The ConvexCombcopula class is used to build convex combinations of copulas,
#' with given positives weights. The rCopula and pCopula functions works for
#' those copulas, assuming they work for the given copulas that we combined
#' in a convex way.
#'
#' See the corresponding vignette for more details about the implementation.
#'
#' @param copulas a list of copulas of same dimension
#' @param alpha a vector of (positive) weights
#' @title Convex Combination of copulas.
#'
#' @name ConvexCombCopula-Class
#' @title Convex Combinations of copulas
#' @rdname ConvexCombCopula-Class
#'
#' @return An instance of the `ConvexCombCopula` S4 class. The object represent the copula that results from a convex combinaison of other copulas, and can be used through several methods to query classical (r/d/p/v)Copula methods, etc.
#' @export
#'
#' @examples
#' dataset <- apply(LifeCycleSavings,2,rank)/(nrow(LifeCycleSavings)+1)
#' copulas <- list(
#'   cbCopula(dataset[,2:3],m=10),
#'   cbCopula(dataset[,2:3],m=5)
#' )
#' alpha <- c(1,4)
#' (cop <- ConvexCombCopula(copulas,alpha))
ConvexCombCopula = function(copulas, alpha = rep(1, length(copulas))) {
  if (missing(copulas) || (!is(copulas, "list"))) {
    if(length(copulas) == 1){
      warning("The copulas argument you provided has only one element. returning this ellement.")
      return(copulas)
    }
    # if(!is(copulas,"Copula")){
    #   stop("The argument copulas must be provided as a list of copulas")
    # } else {
    #
    # }
  }
  if(length(copulas) == 1){
    warning("You provided is a list of only one copula. Returning this copula")
    return(copulas[[1]])
  }
  .ConvexCombCopula(copulas = copulas, alpha = alpha/sum(alpha),dim=copulas[[1]]@dim)
}

setMethod(f = "show",    signature = c(object = "ConvexCombCopula"),                definition = function(object)    {
  cat("This is a ConvexCombCopula , with : \n", "  dim =", dim(object@copulas[[1]]),
      "\n   number of copulas =", length(object@copulas), "\n   alpha =",
      object@alpha, "\n")
  cat("sub-copulas can be accessed trhough the @copulas slot")
})

#' @describeIn rCopula-methods Method for the cbCopula
setMethod(f = "rCopula", signature = c(n = "numeric", copula = "ConvexCombCopula"), definition = function(n, copula) {

  # if n=0, return a 0xdim(copula) matrix :
  if (n == 0) {
    return(matrix(NA, nrow = 0, ncol = dim(copula)))
  }

  # to choose wich copulas will be simulated from, sample
  # 1:length(copulas) with weights equal to alpha, with replacement OFC
  n_cop = length(copula@copulas)
  sampled_copulas <- resample(1:n_cop, size = n, replace = TRUE, prob = copula@alpha)
  tbl <- table(sampled_copulas)

  # then sample from each of those copulas the right number of times :
  samples <- mapply(function(x,y){
    rCopula(n = y, copula = copula@copulas[[x]])
  },as.numeric(names(tbl)),as.vector(tbl),SIMPLIFY = FALSE)

  # then rbind all of them and resample (randomly) rows :
  samples <- do.call(rbind, samples)
  samples <- samples[resample(1:nrow(samples), size = nrow(samples),
                            replace = FALSE), ]
  return(samples)
})

#' @describeIn pCopula-methods Method for the cbCopula
setMethod(f = "pCopula", signature = c(u = "matrix", copula = "ConvexCombCopula"),  definition = function(u, copula) {

  as.vector(
    vapply(copula@copulas,
           pCopula,
           FUN.VALUE=numeric(nrow(u)),
           u=u
    ) %*% copula@alpha
  )

})

