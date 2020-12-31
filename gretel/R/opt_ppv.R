#' Optimize Probabilistic Path Value
#' 
#' Identify the path of optimal probabilistic path value from a source node
#' to a target node.
#' 
#' @param sociomatrix a nonnegative, real valued sociomatrix.
#' @param source an integer index corresponding to a node in \code{sociomatrix}
#' @param target an integer index corresponding to a node in \code{sociomatrix} 
#' @param odds_scale a nonnegative real number indicating the observed tie
#'     strength value that corresponds to 1-1 transmission odds
#' @param odds_scale_by_node sets a transfer odds scale for each node in a 
#'     probabilistic path value calculation.  
#' 
#' @seealso \code{\link{ppv}} to calculate the value of a user-specified path,
#'     \code{\link{all_opt_ppv}} to simultaneously identify the optimal paths 
#'     from any source node to any target node.
#' 
#' @example
#' # Identify the optimal path
#' best_path <- opt_ppv(YangKnoke01, source = 1, target = 4, odds_scale = 1)
#' 
#' # Print the gpv of the optimal path
#' ppv(YangKnoke01, path = best_path, odds_scale = 1)
#' 
#' # Print the gpv of an inferior path
#' ppv(YangKnoke01, path = c(1,2,3,4))
#' 
#' @export
opt_ppv <- function(sociomatrix, source, target, odds_scale = 1, odds_scale_by_node = NULL){
  check_input(sociomatrix, source = source, target = target, odds_scale = odds_scale, odds_scale_by_node = odds_scale_by_node)
  odds_matrix <- scale_to_odds(sociomatrix, odds_scale, odds_scale_by_node)
  probability_matrix <- odds_matrix/(1+odds_matrix)
  lpm <- log(probability_matrix)
  d_eff_mat <- -1*lpm
  path <- shortest_path(d_eff_mat, source, target, node_costs = rep(0,nrow(d_eff_mat)))
  return(path)
}

#' Optimize All Probabilistic Path Values
#' 
#' Identify the path of optimal probabilistic path value from every source 
#' to every target in \code{sociomatrix}.
#' 
#' @param sociomatrix a nonnegative, real valued sociomatrix.
#' @param odds_scale a nonnegative real number indicating the observed tie
#'     strength value that corresponds to 1-1 transmission odds
#' @param odds_scale_by_node sets a transfer odds scale for each node in a 
#'     probabilistic path value calculation.  
#' 
#' @return All optimal paths from source to target nodes in \code{sociomatrix}. To
#'     minimize memory usage, paths are returned as a list of trees in Dijkstra's 
#'     format. Specific paths can be unpacked with \code{unpack} as described in the
#'     example below.
#' 
#' @seealso \code{\link{ppv}} to calculate the value of a user-specified path, 
#'     \code{\link{opt_ppv}} to identify the optimal path from a single source 
#'     node to a single target node
#' 
#' @example 
#' # Identify the optimal paths
#' best_paths <- all_opt_ppv(YangKnoke01, alpha = 1)
#' 
#' # 'best_paths' will contain a list of trees in Dijkstra's format.
#' # 'best_paths[[i]]' is the tree encoding shortest paths from source
#' #     node 'i' to all alters. We can return the optimal path from 
#' #     node 1 to node 4 as follows.
#' unpack(best_paths[[1]], source = 1, target = 4)
#' 
#' @export
all_opt_ppv <- function(sociomatrix, odds_scale = 1, odds_scale_by_node = NULL){
  check_input(sociomatrix, odds_scale = odds_scale, odds_scale_by_node = odds_scale_by_node)
  odds_matrix <- scale_to_odds(sociomatrix, odds_scale, odds_scale_by_node)
  probability_matrix <- odds_matrix/(1+odds_matrix)
  lpm <- log(probability_matrix)
  d_eff_mat <- -1*lpm
  paths <- APSP(d_eff_mat, node_costs = rep(0,nrow(d_eff_mat)))
  # list of spanning trees in dijkstra format
  return(paths)
}