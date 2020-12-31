#' Compute Monte Carlo p-value
#'
#' \code{mc.pvalue} computes the Monte Carlo p-value of each
#' element of \code{tobs} using the relationship
#' \code{(sum(tsim >= x) + 1)/(nsim + 1)} where \code{x} is
#' a specific element of \code{tobs} and \code{nsim =
#' length(tsim)}.
#'
#' @param tobs A vector observed test statistics
#' @param tsim A vector of test statistics from simulated
#'   data
#'
#' @return A vector of p-values.
#' @export
#' @keywords internal
#' @examples
#' mc.pvalue(8:10, 1:9)
mc.pvalue = function(tobs, tsim) {
  nsim = length(tsim)
  unname(sapply(tobs, function(x) (sum(tsim >= x) + 1) / (nsim + 1)))
}
