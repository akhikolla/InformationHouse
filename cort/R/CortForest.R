.CortForest = setClass(Class = "CortForest", contains = "empiricalCopula",
                 slots = c(p_value_for_dim_red = "numeric",
                           number_max_dim = "numeric",
                           verbose_lvl="numeric",
                           trees = "list",
                           weights = "numeric",
                           indexes = "matrix",
                           pmf = "matrix",
                           norm_matrix = "matrix",
                           oob_pmf = "matrix",
                           oob_kl = "numeric",
                           oob_ise = "numeric"),
                 validity = function(object) {
                   errors <- c()
                   if(!(length(dim(object@data))==2)){
                     errors <- c(errors, "data should be a matrix")
                   }
                   if(!(all(object@data >= 0) & all(1 > object@data))){
                     errors <- c(errors,"data should be a matrix inside 0,1")
                   }
                   if(object@number_max_dim > ncol(object@data)){
                     errors <- c(errors, "the maximum number of splitting dimensions should be smaller than the number of dimensons!")
                   }
                   if (length(errors) == 0)
                     TRUE else errors
                 })

#' CortForest class
#'
#' This class implements the bagging of CORT models, with an out-of-bag error minimisation in the weights.
#'
#'
#' See O. Laverny, V. Maume-Deschamps, E. Masiello and D. Rullière (2020) for the details of this density estimation procedure, and `vignettes(package='cort')` for examples of usecases.
#'
#'
#' @param x The data, must be provided as a matrix with each row as an observation.
#' @param p_value_for_dim_red a p_value for the localised dimension reduction test
#' @param min_node_size The minimum number of observation avaliable in a leaf to initialise a split.
#' @param pseudo_data set to True if you are already providing data on the copula space.
#' @param compte_loo_weights Defaults to FALSE. Allows to use an automatic re-weighting of the trees in the forest, based on leave-one-out considerations.
#' @param number_max_dim The maximum number of dimension a split occurs in. Defaults to be all of the dimensions.
#' @param n_trees Number of trees
#' @param verbose_lvl verbosity level : can be 0 (default) or an integer. bigger the integer bigger the output level.
#' @param force_grid boolean (default: FALSE). set to TRUE to force breakpoint to be on the n-checkerboard grid in every tree.
#' @param oob_weighting boolean (default : TRUE) option to weight the trees with an oob criterion (otherwise they are equally weighted)
#'
#' @name CortForest-Class
#' @title Bagged Cort copulas
#' @rdname CortForest-Class
#'
#' @return An instance of the `CortForest` S4 class. The object represent the fitted copula and can be used through several methods to query classical (r/d/p/v)Copula methods, constraint influence, etc.
#' Beside returning some inputted parameters, notable slots are :
#'
#' \itemize{
#'  \item{`trees` }{A list of Cort objects representing each fitted tree in the forest.}
#'  \item{`weights` }{The weigths of each tree.}
#'  \item{`indexes` }{The indexes of data points that were selected for fitting the trees}
#'  \item{`pmf` }{The density of each tree on data points}
#'  \item{`norm_matrix` }{The matrix of scalar product between trees}
#'  \item{`oob_pmf` }{The density of each tree on data points it did not see during fitting}
#'  \item{`oob_kl` }{The out-of-bag Kullback-Leibler divergence of each tree }
#'  \item{`oob_ise` }{The out-of-bag Integrated Square Error of each tree}
#' }
#'
#' More details about these slots can be found in the reference.
#'
#' @export
#'
#' @references
#' \insertRef{laverny2020}{cort}
#'
#' @examples
#' (CortForest(LifeCycleSavings[,1:3],number_max_dim=2,n_trees=2))
CortForest = function(x,
                p_value_for_dim_red=0.75,
                n_trees = 10,
                compte_loo_weights = FALSE,
                min_node_size=1,
                pseudo_data=FALSE,
                number_max_dim=NULL,
                verbose_lvl=2,
                force_grid = FALSE,
                oob_weighting = TRUE) {

  # Furrr options to set seed :
  FO = furrr::furrr_options(seed=TRUE)

  # Deal with verbosity :
  showing = ifelse(verbose_lvl <= 1,"","========================")

  # Extract and coerce data :
  data= as.matrix(x)
  row.names(data) = NULL
  colnames(data) = NULL
  if(!pseudo_data){ data = apply(data,2,rank,ties.method="random")/(nrow(data)+1) }
  d = ncol(data)
  n = nrow(data)

  # Bootstrap observation for the trees
  indexes = matrix(resample(1:n,size=n*n_trees,replace = TRUE),nrow=n_trees,ncol=n) # n_trees, n


  # Fit the trees
  if(verbose_lvl>0){cat(showing,"Computing trees...\n")}
  trees = furrr::future_map(1:n_trees,function(i){
    Cort(data[indexes[i,],],
         p_value_for_dim_red=p_value_for_dim_red,
         min_node_size=min_node_size,
         pseudo_data=TRUE,
         number_max_dim=number_max_dim,
         verbose_lvl=0,
         force_grid = force_grid)},.progress=TRUE,.options = FO)

  # Compute the stats
  if(verbose_lvl>0){cat(showing,"Computing statistics...\n")}

  # Computation of the pmf and of it's mask
  if(verbose_lvl>1){cat("     Computing pmf...\n")}
  pmf = is_out = matrix(NA,nrow=n_trees,ncol=n)
  for (i in 1:n_trees){
    is_out[i,] = 1-(1:n %in% indexes[i,])
    pmf[i,] = dCopula(data,trees[[i]])
  }

  # Parallel computation for norm_matrix :
  if(verbose_lvl>1){cat("     Computing norm matrix...\n")}
  norm_matrix = diag(purrr::map_dbl(trees,quad_norm))/2
  idx = data.frame(nrow = rep(1:n_trees,n_trees),ncol=rep(1:n_trees,each=n_trees))
  idx = idx[idx$nrow < idx$ncol,]
  norm_matrix[as.matrix(idx)] <- furrr::future_map_dbl(1:nrow(idx),function(i){quad_prod(trees[[idx[i,1]]],trees[[idx[i,2]]])},.progress=TRUE,.options = FO)
  norm_matrix = norm_matrix + t(norm_matrix)

  # Fitting weights
  if(!oob_weighting){
    oob_wts = rep(1/n_trees,n_trees)
  } else {
    if(verbose_lvl>1){cat("     Computing weights...\n")}

      loss <- function(opt_par, pmf, norm_mat, is_out){
        weights_out = is_out*opt_par
        return(opt_par %*% norm_mat %*% opt_par - mean(colSums(pmf*weights_out)/colSums(weights_out),na.rm=TRUE))
      }

      # a good starting point :
      x0 = 1/diag(norm_matrix)
      x0 = x0 / sum(x0)

      rez = nloptr::slsqp(
        x0 = x0,
        fn = loss,
        lower = rep(0,n_trees),
        heq = function(opt_par){sum(opt_par)-1},
        nl.info=verbose_lvl>1,
        control=list(maxeval  = 1000000L,
                     xtol_rel = 1e-10,
                     ftol_rel = 1e-10,
                     ftol_abs = 1e-10),
        pmf = pmf,
        norm_mat = norm_matrix,
        is_out = is_out)$par

      oob_wts = pmax(rez,0)/sum(pmax(rez,0))
  }

  # OComputation of out-of-bag statistics
  if(verbose_lvl>1){cat("     Computing oob stats...\n")}
  oob_pmf = matrix(0,nrow=n,ncol=n_trees - 1)
  oob_kl = numeric(n_trees - 1)
  oob_ise = numeric(n_trees - 1)
  oob_wts_out = is_out*oob_wts
  for (j in 2:n_trees){
    oob_pmf[,j-1] = colSums(pmf[1:j,]*oob_wts_out[1:j,])/colSums(oob_wts_out[1:j,])
    oob_kl[j-1] = -mean(log(oob_pmf[,j-1]),na.rm=TRUE)
    oob_ise[j-1] = (oob_wts[1:j] %*% norm_matrix[1:j,1:j] %*% oob_wts[1:j])/(sum(oob_wts[1:j])^2) -2 * mean(oob_pmf[,j-1],na.rm=TRUE)
  }

  if(verbose_lvl>0){cat(showing,"Done !\n")}
  return(.CortForest(
    data = data,
    p_value_for_dim_red = p_value_for_dim_red,
    number_max_dim = max(number_max_dim,d),
    verbose_lvl=verbose_lvl,
    trees = trees,
    dim = ncol(data),
    weights = oob_wts,
    indexes = indexes,
    pmf = pmf,
    norm_matrix = norm_matrix,
    oob_pmf = oob_pmf,
    oob_kl = oob_kl,
    oob_ise = oob_ise
  ))
}

setMethod(f = "show", signature = c(object = "CortForest"), definition = function(object){
  cat(paste0("CortForest copula model: ",nrow(object@data),"x",ncol(object@data),"-dataset and ",
             length(object@trees)," Cort trees."))
})

#' @describeIn rCopula-methods Method for the class CortForest
setMethod(f = "rCopula", signature = c(n = "numeric", copula = "CortForest"), definition = function(n, copula) {
  if(n==0){
    return(matrix(0,nrow=0,ncol=copula@dim))
  }
  sampled_indexes = resample(1:length(copula@trees),size=n,prob = copula@weights,replace=TRUE)
  return(do.call(rbind,purrr::map(unique(sampled_indexes),~rCopula(sum(sampled_indexes==.x),copula@trees[[.x]]))))

})

#' @describeIn pCopula-methods Method for the class CortForest
setMethod(f = "pCopula", signature = c(u = "matrix",  copula = "CortForest"), definition = function(u, copula) {
  return(as.vector(vapply(copula@trees,function(t){pCopula(u,t)},rep(0.5,nrow(u))) %*% copula@weights))
})

#' @describeIn dCopula-methods Method for the class CortForest
setMethod(f = "dCopula", signature = c(u = "matrix",  copula="CortForest"),   definition = function(u, copula) {
  return(vapply(copula@trees,function(t){dCopula(u,t)},rep(0.5,nrow(u))) %*% copula@weights)
})

#' @describeIn constraint_infl-methods Method for the class CortForest
setMethod(f = "constraint_infl", signature = c(object="CortForest"),   definition = function(object) {
  return(purrr::map_dbl(object@trees,constraint_infl))
})

#' @describeIn quad_norm-methods Method for the class CortForest
setMethod(f = "quad_norm", signature = c(object="CortForest"),   definition = function(object) {
  return(object@weights %*% object@norm_matrix %*% object@weights)
})













