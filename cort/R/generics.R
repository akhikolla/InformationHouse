#' @import methods
#' @importFrom Rdpack reprompt
#' @importFrom stats runif
#' @importFrom utils capture.output
#' @include utils.R
NULL

#' Copula volume on hyper-boxes
#'
#' u must be piecewise smaller than v, otherwise the function will return an error.
#'
#' A method is currently implemented for the main virtual class 'Copula', but it assumes
#' that a pCopula method is avaliable for the given copula. This method could be used with Copulas that are not from this package, assuming that pCopula(u,cop) works.
#'
#' This function computes the measure of the copula according to the algorithm proposed by Cherubini U, Romagnoli S (2009-oct).
#'
#'
#' @param u numeric matrix : minimum point of the hyper-rectangles, one row per observation.
#' @param v numeric matrix : maximum point of the hyper-rectangle, one row per observation.
#' @param copula the copula that we compute the measure on the box (u,v)
#' @param ... other parameter to be passed to methods for this generic.
#'
#' @return the measure of the copula.
#' @exportMethod vCopula
#' @name vCopula
#' @rdname vCopula-methods
#'
#' @examples
#' cop <- cbCopula(LifeCycleSavings,m = 5)
#' vCopula(rep(0,5),rep(1,5),cop) == 1
#' vCopula(rep(0,5),rep(0.5,5),cop)
#'
#' @references
#' \insertRef{cherubini2009}{cort}
#'
setGeneric("vCopula", function(u, v, copula, ...) {

  # taken from the generic of pCopula, does mainly the same...

  u <- normalise_data(u,copula@dim)
  v <- normalise_data(v,copula@dim)
  standardGeneric("vCopula")
})


#' Copula density
#'
#' This function returns the density of a given copula on given observations.
#'
#'
#' @param u numeric matrix : one row per observation
#' @param copula the copula object
#' @param ... other parameter to be passed to methods for this generic.
#'
#' @return The density of the copula on each observation
#' @exportMethod dCopula
#' @name dCopula
#' @rdname dCopula-methods
#'
#' @examples
#' cop <- cbCopula(cort::funcdep_data[1:10,1:2], m = 5)
#' dCopula(rep(0,2),cop)
#' dCopula(rep(0.5,2),cop)
#' dCopula(rep(1,2),cop)
#'
setGeneric("dCopula", function(u, copula, ...) {
  u <- normalise_data(u,copula@dim)
  standardGeneric("dCopula")
})

#' Copula cdf
#'
#' This function returns the value of the copula itself on given points.
#'
#'
#' @param u numeric matrix : one row per observation
#' @param copula the copula object
#' @param ... other parameter to be passed to methods for this generic.
#'
#' @return The value of the copula on each observation
#' @exportMethod pCopula
#' @name pCopula
#' @rdname pCopula-methods
#'
#' @examples
#' cop <- cbCopula(cort::recoveryourself_data,m = 5)
#' pCopula(rep(0,2),cop) == 0
#' pCopula(rep(0.5,2),cop)
#' pCopula(rep(1,2),cop) == 1
#'
setGeneric("pCopula", function(u, copula, ...) {
  u <- normalise_data(u,copula@dim)
  standardGeneric("pCopula")
})

#' @rdname vCopula-methods
#' @aliases vCopula,matrix,matrix,Copula
setMethod("vCopula", signature = c(u = "matrix", v = "matrix"),
          definition = function(u, v, copula) {

            # can handle any copula thant pCopula could handle.
            # shoul be better vetorised...
            # u and v must be numeric, copula must be a copula, and v must be
            # smaller than u
            if (nrow(u) != nrow(v)) {
              stop("u and v must have same shape (same number of row and columns)")
            }
            if((nrow(u) == 0) || (nrow(v) == 0)){
              return(numeric())
            }

            if (any(v < u)) {
              stop("u must be smaller than v !")
            }

            # fastening ? maybe not...
            pCop_method <- selectMethod(pCopula, c("matrix", class(copula)))


            d = dim(copula)
            p <- t(sapply(1:(2^d),function(i){number2binary(i-1,d)}))
            sign <- (-1)^rowSums(p)
            return(sapply(1:nrow(u),function(i){
              if(all(u[i,] == v[i,])){return(0)}
              eval_points <-t(t(p) * as.vector(u[i,]) + t(1-p) * as.vector(v[i,]))
              return(sum(sign * pCopula(eval_points,copula)))
            }))
          })


#' Copula random generation
#'
#' Random number generation following the given copula. This function performs the simulation of random vectors following the copula.
#'
#'
#' @param n the number of simulations
#' @param copula the copula object
#' @param ... other parameter to be passed to methods for this generic.
#'
#' @return A matrix with `n` rows, each representing a random vector generated from the provided copula.
#' @exportMethod rCopula
#' @name rCopula
#' @rdname rCopula-methods
#'
#' @examples
#' cop <- cbCopula(cort::clayton_data,m = 5)
#' xx <- rCopula(1000,cop)
#'
setGeneric("rCopula", function(n, copula, ...) standardGeneric("rCopula"))


#' Spearman's rho matrix of a copula
#'
#' Computes the bivariate Spearmann's rho matrix for a copula.
#'
#'
#' @param copula the copula object
#'
#' @return the density of the copula on each observation
#' @exportMethod biv_rho
#' @name biv_rho
#' @rdname biv_rho-methods
#'
#' @examples
#' cop <- Cort(LifeCycleSavings[,1:3])
#' biv_rho(cop)
#'
setGeneric("biv_rho", function(copula) standardGeneric("biv_rho"))

#' Kendall's tau matrix of a copula
#'
#' Computes the bivariate Kendall's tau matrix for a copula.
#'
#'
#' @param copula the copula object
#'
#' @return the density of the copula on each observation
#' @exportMethod biv_tau
#' @name biv_tau
#' @rdname biv_tau-methods
#'
#' @examples
#' cop <- Cort(cort::funcdep_data[1:10,1:3])
#' biv_tau(cop)
#'
setGeneric("biv_tau", function(copula) standardGeneric("biv_tau"))




#' Loss of a copula estimation (if the model has one)
#'
#' Currently only implemented for Cort models.
#' Compute the loss of the model.
#'
#'
#' @param object the copula object
#'
#' @return the Integrated square error loss of the model
#' @exportMethod loss
#' @name loss
#' @rdname loss-methods
#'
#' @examples
#' cop <- Cort(cort::recoveryourself_data[1:10,])
#' loss(cop)
#'
setGeneric("loss", function(object) standardGeneric("loss"))


#' Constraint influence of the model (if it has one)
#'
#' Currently only implemented for Cort models.
#' Compute the constraint influence of the model
#'
#'
#' @param object the copula object
#'
#' @return The constraint influence statistic of the model
#' @exportMethod constraint_infl
#' @name constraint_infl
#' @rdname constraint_infl-methods
#'
#' @examples
#' cop <- Cort(cort::recoveryourself_data[1:10,])
#' constraint_infl(cop)
#'
setGeneric("constraint_infl", function(object) standardGeneric("constraint_infl"))

#' Quadratic norm of the model (if it has one)
#'
#' Currently only implemented for Cort models.
#' Compute the L2 norm of the model
#'
#'
#' @param object the copula object
#'
#' @return the Integrated square error quad_norm of the model
#' @exportMethod quad_norm
#' @name quad_norm
#' @rdname quad_norm-methods
#'
#' @examples
#' cop <- Cort(cort::impossible_data)
#' quad_norm(cop)
#'
setGeneric("quad_norm", function(object) standardGeneric("quad_norm"))

#' Quadratic product with data of the model (if it has one)
#'
#' Currently only implemented for Cort models.
#' Compute the quadratic product with the empirical density from the data
#'
#'
#' @param object the copula object
#'
#' @return the quad_prod_with_data of the model
#' @exportMethod quad_prod_with_data
#' @name quad_prod_with_data
#' @rdname quad_prod_with_data-methods
#'
#' @examples
#' cop <- Cort(LifeCycleSavings[,1:3])
#' quad_prod_with_data(cop)
#'
setGeneric("quad_prod_with_data", function(object) standardGeneric("quad_prod_with_data"))

#' Quadratic product of two copulas (if they have one)
#'
#' Currently only implemented for Cort models.
#' Compute the L2 quadratic product of 2 trees
#'
#'
#' @param object : the tree
#' @param other_tree : the other tree
#'
#' @return the quadratic product between the trees
#' @exportMethod quad_prod
#' @name quad_prod
#' @rdname quad_prod-methods
#'
#' @examples
#' cop <- Cort(LifeCycleSavings[,1:3])
#' all.equal(quad_prod(cop,cop),quad_norm(cop))
#'
setGeneric("quad_prod", function(object,other_tree) standardGeneric("quad_prod"))


#' Kendall function of a copula (if it has one)
#'
#' Currently only implemented for Cort models.
#' Compute the Kendall cdf from the model in a point t
#'
#'
#' @param object : the tree
#' @param t : the value where to compute the kendall function, may be a vector of evaluation values;
#' @param ... other parameters passed to methods
#'
#' @return the quadratic product between the trees
#' @exportMethod kendall_func
#' @name kendall_func
#' @rdname kendall_func-methods
#'
#' @examples
#' cop <- Cort(LifeCycleSavings[,1:3])
#' kendall_func(cop,0.5)
#'
setGeneric("kendall_func", function(object,t,...) {
  standardGeneric("kendall_func")
})

#' Projection on smaller dimensions of a copula (if implemented)
#'
#' Currently only implemented for Cort models.
#' Compute, as a Cort object, the projection on a smaller set of dimensions of a Cort object.
#'
#'
#' @param object : the tree
#' @param dims the set of dimensions
#'
#' @return other cort object
#' @exportMethod project_on_dims
#' @name project_on_dims
#' @rdname project_on_dims-methods
#'
#' @examples
#' cop <- Cort(LifeCycleSavings[,1:3])
#' projection = project_on_dims(cop,c(1,2))
#'
setGeneric("project_on_dims", function(object,dims) standardGeneric("project_on_dims"))


