#' Calculate probabilistic path value
#'
#' Given a real valued sociomatrix, a path, and an optional \code{odds_scale}, \code{ppv}
#' calculates the transmission odds for the path and returns the transmission odds
#' times \code{odds_scale} so the result can be directly compared with observed 
#' tie strenghts.
#'
#' We assume that observed tie strengths in \code{sociomatrix} are linearly
#' proportional to transmission odds. That is, if the transmission odds for
#' a strength 1 tie are 1 to 1, the transmission odds for a strength 5 tie are
#' 1 to 5.
#'
#'
#' @param sociomatrix a nonnegative, real valued sociomatrix.
#' @param path an integer vector of node indices from \code{sociomatrix}.
#' @param odds_scale a nonnegative real number indicating the observed tie
#'     strength value that corresponds to 1-1 transmission odds
#' @param odds_scale_by_node sets a transfer odds scale for each node in a 
#'     probabilistic path value calculation.  
#'     
#' @seealso \code{\link{opt_ppv}} to identify the path of optimal 'ppv' between two nodes 
#'     and \code{\link{all_opt_ppv}} to identify the optimal paths between all pairs of 
#'     nodes. Calling \code{\link{generate_proximities}} with \code{mode = 'ppv'} 
#'     returns a matrix 'ppv' values for the optimal paths between all pairs of 
#'     nodes.
#' 
#' @examples
#' ## Calculate ppv along a path in a sociomatrix
#' ppv(YangKnoke01, path = c(1,2,5), odds_scale = 3)
#' 
#' ## This path doesn't exist
#' gpv(YangKnoke01, path = c(1,2,4,5))
#' 
#' @export
ppv <- function(sociomatrix, path, odds_scale = 1, odds_scale_by_node = NULL){
  check_input(sociomatrix, path = path, odds_scale = odds_scale, odds_scale_by_node = odds_scale_by_node)
  odds_matrix <- scale_to_odds(sociomatrix, odds_scale, odds_scale_by_node)
  probability_matrix <- odds_matrix/(1+odds_matrix)
  lpm <- log(probability_matrix)
  log_prob <- sum_path(lpm, path)
  transmission_prob <- exp(log_prob)
  transmission_odds <- transmission_prob/(1 - transmission_prob)
  return(transmission_odds*odds_scale)
}