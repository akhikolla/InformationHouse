#' Pairwise Tests for Association/Correlation Between Paired Samples
#'
#' Calculate pairwise correlations between group levels with corrections for multiple testing.
#'
#' Note that correlation tests require that the two vectors examined are of the same length.
#' Thus, if the grouping defines groups of varying lengths then the specific correlation is
#' not computed and a \code{NA} is returned instead. The adjusted p values are only based on
#' the actual correlation that are computed.
#' Extra arguments that are passed on to \code{cor.test} may or may not be sensible in this context.
#' 
#' @param x response vector.
#' @param g grouping vector or factor.
#' @param p.adjust.method method for adjusting p values (see \code{\link[stats]{p.adjust}}). Can be abbreviated.
#' @param method string argument to set the method to compute the correlation. Possibilities are "pearson" (the default), "kendall", and "spearman"
#' @param ... additional arguments passed to \code{\link[stats]{cor.test}}.
#'
#' @return Object of class \code{pairwise.htest}
#'
#' @examples
#'
#' attach(airquality)
#' Month <- factor(Month, labels = month.abb[5:9])
#' pairwise.cor.test(Ozone, Month)
#' pairwise.cor.test(Ozone, Month, p.adj = "bonf")
#' detach()
#' @export
pairwise.cor.test <- function (x, g, p.adjust.method = p.adjust.methods,
                               method = c("pearson", "kendall", "spearman"), ...)
{
  method <- match.arg(method)
  p.adjust.method <- match.arg(p.adjust.method)
  DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(g)))
  g <- factor(g)
  compare.levels <- function(i, j) {
    xi <- x[as.integer(g) == i]
    xj <- x[as.integer(g) == j]

    # The correlations require that the groups are of the same length
    if (length(xi) != length(xj)) {
        warning("Group sizes are not identical")
        return(NA)
    } else {
        cor.test(xi, xj, method=method, ...)$p.value
    }
  }
  PVAL <- stats::pairwise.table(compare.levels, levels(g), p.adjust.method)
  if (method=="pearson")
    METHOD <- "Pearson's product-moment correlation"
  if (method=="kendall")
    METHOD <- "Kendall's rank correlation tau"
  if (method=="spearman")
    METHOD <- "Spearman's rank correlation rho"

  ans <- list(method = METHOD, data.name = DNAME, p.value = PVAL,
              p.adjust.method = p.adjust.method)
  class(ans) <- "pairwise.htest"
  ans
}
