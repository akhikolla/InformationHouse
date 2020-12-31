#' Return most significant, non-overlapping zones
#'
#' \code{sig_noc} return the significant, non-overlapping
#' zones order from most significant to least significant.
#' @param order_by Either \code{'tobs'} or \code{'pvalue'},
#' indicating the argument by which to order \code{zones}.
#' @inheritParams smerc_cluster
#'
#' @return A list with the significant, ordered,
#' non-overlapping \code{tobs}, \code{zones}, \code{pvalue}.,
#' and \code{idx} (a vector with the relevant indices of
#' the original zones).
#' @export
#' @examples
#' tobs = c(1, 3, 2)
#' zones = list(1:2, 1:3, 2:3)
#' pvalue = c(0.5, 0.01, 0.02)
#' sig_noc(tobs, zones, pvalue, alpha = 0.05)
sig_noc = function(tobs, zones, pvalue, alpha,
                   order_by = "tobs") {
  # argument checking
  N = length(tobs)
  arg_check_tobs(tobs)
  arg_check_zones(zones, N)
  arg_check_pvalue(pvalue, N)
  arg_check_alpha(alpha)
  if (length(order_by) != 1) {
    stop("order_by must have length 1")
  }
  if (!is.element(order_by, c("tobs", "pvalue"))) {
    stop("order_by must be 'tobs' or 'pvalue'")
  }

  # create idx from original list of zones
  idx = seq_len(N)
  # determine if there are any significant zones
  minp = which.min(pvalue)
  if (pvalue[minp] > alpha) {
    warning("No significant clusters.  Returning most likely cluster.")
    sig = minp
  } else {
    sig = which(pvalue <= alpha)
  }
  # only keep significant zones and info
  tobs = tobs[sig]
  zones = zones[sig]
  pvalue = pvalue[sig]
  idx = idx[sig]

  # order zones from most to least significant
  if (order_by == "tobs") {
    ozones = order(tobs, decreasing = TRUE)
  } else {
    ozones = order(pvalue, decreasing = FALSE)
  }
  # reorder zones and info
  tobs = tobs[ozones]
  zones = zones[ozones]
  pvalue = pvalue[ozones]
  idx = idx[ozones]

  # determine significant non-overlapping clusters
  # in order of significance
  # sig = noz(zones)
  # use c++ for substantial efficiency increase
  # offset by 1 since c++ starts index at 0
  sig = noc_cpp(zones) + 1
  return(list(tobs = tobs[sig],
              zones = zones[sig],
              pvalue = pvalue[sig],
              idx = idx[sig]))
}






