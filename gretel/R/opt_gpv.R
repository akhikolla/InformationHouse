#' Optimize Generalized Path Value
#' 
#' Identify the path of optimal generalized path value from a source node
#' to a target node.
#' 
#' @param sociomatrix a nonnegative, real valued sociomatrix.
#' @param source an integer index corresponding to a node in \code{sociomatrix}
#' @param target an integer index corresponding to a node in \code{sociomatrix}
#' @param p a nonnegative real number that sets the 'p-norm' parameter for 
#'     generalized path value calculation.
#' @param node_costs a list of costs, in order, of all nodes represented in the 
#'     sociomatrix, all are assumed 0 if unspecified
#'
#' @seealso \code{\link{gpv}} to calculate the value of a user-specified path,
#'     \code{\link{all_opt_gpv}} to simultaneously identify the optimal paths 
#'     from any source node to any target node.
#' 
#' @example
#' # Identify the optimal path
#' best_path <- opt_gpv(YangKnoke01, source = 1, target = 4, p = 1)
#' 
#' # Print the gpv of the optimal path
#' gpv(YangKnoke01, path = best_path, p = 1)
#' 
#' # Print the gpv of an inferior path
#' gpv(YangKnoke01, path = c(1,2,3,4))
#' 
#' @export
opt_gpv <- function(sociomatrix, source, target, p = Inf, node_costs = NULL){
  check_input(sociomatrix, source = source, target = target, p_norm = p, node_costs = node_costs)
  if(is.null(node_costs)){
    node_costs <- rep(0,nrow(sociomatrix))
  }
  if(p == 0){
    d_eff_mat <- 1/as_binary(sociomatrix)
    path <- shortest_path(d_eff_mat, source, target, node_costs)
  }else if(p < Inf){
    d_eff_mat <- 1/(sociomatrix^p)
    path <- shortest_path(d_eff_mat, source, target, node_costs)
  }else if(p == Inf){
    d_eff_mat <- 1/sociomatrix
    path <- shortest_path(d_eff_mat, source, target, node_costs = NULL, p_finite = F)
  }
  return(path)
}

#' Optimize All Generalized Path Values
#' 
#' Identify the path of optimal generalized path value from every source 
#' to every target in \code{sociomatrix}.
#' 
#' @param sociomatrix a nonnegative, real valued sociomatrix.
#' @param p a nonnegative real number that sets the 'p-norm' parameter for 
#'     generalized path value calculation.
#' @param node_costs a list of costs, in order, of all nodes represented in the 
#'     sociomatrix, all are assumed 0 if unspecified
#' 
#' @return All optimal paths from source to target nodes in \code{sociomatrix}. To
#'     minimize memory usage, paths are returned as a list of trees in Dijkstra's 
#'     format. Specific paths can be unpacked with \code{unpack} as described in the
#'     example below.
#' 
#' @seealso \code{\link{gpv}} to calculate the value of a user-specified path, 
#'     \code{\link{opt_gpv}} to identify the optimal path from a single source 
#'     node to a single target node
#' 
#' @example 
#' # Identify the optimal paths
#' best_paths <- all_opt_gpv(YangKnoke01, p = 1)
#' 
#' # 'best_paths' will contain a list of trees in Dijkstra's format.
#' # 'best_paths[[i]]' is the tree encoding shortest paths from source
#' #     node 'i' to all alters. We can return the optimal path from 
#' #     node 1 to node 4 as follows.
#' unpack(best_paths[[1]], source = 1, target = 4)
#' 
#' @export
all_opt_gpv <- function(sociomatrix, p = Inf, node_costs = NULL){
  check_input(sociomatrix, p_norm = p, node_costs = node_costs)
  if(is.null(node_costs)){
    node_costs <- rep(0,nrow(sociomatrix))
  }
  if(p == 0){
    d_eff_mat <- 1/as_binary(sociomatrix)
    paths <- APSP(d_eff_mat, node_costs)
  }else if(p < Inf){
    d_eff_mat <- 1/(sociomatrix^p)
    paths <- APSP(d_eff_mat, node_costs)
  }else if(p == Inf){
    d_eff_mat <- 1/sociomatrix
    paths <- APSP(d_eff_mat, node_costs = NULL, p_finite = F)
  }
  # list of spanning trees in dijkstra format
  return(paths)
}