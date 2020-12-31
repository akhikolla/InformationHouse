#' Generate a Proximity Matrix
#' 
#' Generates a proximity matrix in one of three modes:
#' \describe{
#'   \item{\code{'ogpv'}}{Optimal Generalized Path Value. Entry \code{i,j} of the 
#'     proximity matrix will equal the optimal 'gpv' among all paths connecting node 
#'     \code{i} to node \code{j}.}
#'   \item{\code{'oppv'}}{Optimal Probabilistic Path Value. Entry \code{i,j} of the 
#'     proximity matrix will equal the optimal 'ppv' among all paths connecting node 
#'     \code{i} to node \code{j}.}
#'   \item{\code{'sconductivity'}}{Social Conductivity (Random Walk Probability). If 
#'     each tie strength recorded in \code{sociomatrix} is taken to be analogous to the 
#'     conductivity of an electrical component, \code{i,j} of the proximity matrix 
#'     will equal total conductivity of all paths from node \code{i} to node \code{j}.}
#' }
#' 
#' @param sociomatrix a nonnegative, real valued sociomatrix.
#' @param mode a selection of \code{'ogpv'}, \code{'oppv'}, or \code{'sconductivity'}
#' @param p if \code{mode} is \code{'ogpv'}, determines 'p-norm' parameter for 
#'     generalized path value calculation.
#' @param node_costs if \code{mode} is \code{'ogpv'}, assigns transmission costs to
#'     vertices within the sociomatrix.
#' @param odds_scale if \code{mode} is \code{'oppv'}, sets a global transfer odds 
#'     scale for probabilistic path value calculation.
#' @param odds_scale_by_node if \code{mode} is \code{'oppv'}, sets a transfer odds 
#'     scale for each node in a probabilistic path value calculation.
#' 
#' @examples
#' ## Generate a proximity matrix in each mode
#' ## Optimal Generalized Path Value
#' generate_proximities(YangKnoke01, mode = "ogpv", p = Inf, node_costs = c(1,3,3,2,1))
#' 
#' ## Optimal Probabilistic Path Value
#' generate_proximities(YangKnoke01, mode = "oppv", odds_scale = 2)
#' 
#' ## Sconductivity
#' generate_proximities(YangKnoke01, mode = "sconductivity")
#' 
#' @seealso \code{\link{gpv}}, \code{\link{ppv}}
#' @export
generate_proximities <- function(sociomatrix, mode = c("ogpv", "oppv", "sconductivity"), p = Inf, node_costs = NULL, odds_scale = 1, odds_scale_by_node = NULL){
  check_input(sociomatrix, p_norm = p, node_costs = node_costs, odds_scale = odds_scale, odds_scale_by_node = odds_scale_by_node)
  mode <- mode[1] # "ogpv" is default if no mode is specified.
  if(!(mode %in% c("ogpv", "oppv", "sconductivity"))){
    stop("'mode' must be one of 'ogpv', 'oppv', or 'sconductivity'")
  }
  nv <- nrow(sociomatrix)
  if(mode == "ogpv"){
    if((!is.null(odds_scale_by_node)) || (odds_scale != 1)){warning("Odds scale parameters not used for 'mode = ogpv'")}
    trees <- all_opt_gpv(sociomatrix, p, node_costs)
    proximity_matrix <- matrix(0, nrow = nv, ncol = nv)
    for(i in 1:nv){
      tree <- trees[[i]]
      for(j in 1:nv){
        path <- unpack(tree,i,j)
        if(i == j){
          proximity_matrix[i,j] <- Inf
        }else if(anyNA(path)){
          proximity_matrix[i,j] <- 0
        }else{
          proximity_matrix[i,j] <- gpv(sociomatrix, path, p, node_costs)
        }
      }
    }
  }else if(mode == "oppv"){
    if((!is.null(node_costs)) || (p != Inf)){warning("Node cost parameter and p norm parameter not used for 'mode = oppv'")}
    trees <- all_opt_ppv(sociomatrix, odds_scale, odds_scale_by_node)
    proximity_matrix <- matrix(0, nrow = nv, ncol = nv)
    for(i in 1:nv){
      tree <- trees[[i]]
      for(j in 1:nv){
        path <- unpack(tree,i,j)
        if(i == j){
          proximity_matrix[i,j] <- Inf
        }else if(anyNA(path)){
          proximity_matrix[i,j] <- 0
        }else{
          proximity_matrix[i,j] <- ppv(sociomatrix, path, odds_scale, odds_scale_by_node)
        }
      }
    }
  }else if(mode == "sconductivity"){
    if((!is.null(node_costs)) || (p != Inf) || (!is.null(odds_scale_by_node)) || (odds_scale != 1)){
      warning("Odds scale, node cost, and p norm parameters not used for 'mode = sconductivity'")
    }
    if(!isSymmetric(sociomatrix)){
      warning("sconductivity beta - intended only for undirected networks and no node costs")
    }
    proximity_matrix <- sconduct(sociomatrix)
  }
  return(proximity_matrix)
}