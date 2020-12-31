#' Plot the outcome of a paired equivalence test
#'
#' @param x the object to be plotted
#' @param shade logical. Should the results be plotted using a shaded
#'   equivalence region?
#' @param ... arguments passed to \code{ggplot2::theme}.
#' @param results data frame. The \code{results} component of a
#'   \code{paired_equivalence} object
#'
#' @return A plot of the equivalence test
#' @export
#'
#' @details
#'
#' \code{shaded_equivalence_plot} plots the results of an equivalence test in
#' which a single equivalence region applies to all variables. In that case, the
#' equivalence region is displayed as a shaded region.
#' \code{unshaded_equivalence_plot} plots the results of an equivalence test in
#' which variables have unique equivalence regions. In that case, the
#' equivalence regions are displayed as dodged "confidence intervals".
#'
#'
#' @examples
#' set.seed(1544)
#' y <- rnorm(500, 17.4, 9)
#' z <- data.frame(
#'   var1 = rnorm(500, 15, 4),
#'   var2 = rnorm(500, 23, 7.3)
#' )
#'
#' # Optionally create artificial missing values to trigger unshaded plot
#' missing_indices <- sample(seq(nrow(z)), 250)
#' z$var1[missing_indices] <- NA
#'
#' x <- paired_equivalence_test(
#'   z, y, "criterion", scale = "relative",
#'   relative_region_width = 0.25
#' )
#'
#' plot(x)
plot.paired_equivalence <- function(x, shade = "auto", ...) {

  results <- x$results
  exclude <- results$y_type == "comparison" &
    results$scale == "relative"

  stopifnot(!all(exclude))
  if (any(exclude)) warning(paste(
    "CI plotting not possible when",
    "y_type == \"comparison\" and\n  ",
    "scale == \"relative\".",
    "Removing", sum(exclude), "results."
  ))

  results <- results[!exclude, ]

  if (shade == "auto") {

    if (nrow(results) == 1) {
      shade <- TRUE
    } else {
      shade <- all(
        duplicated(
          results[ ,c("region_low", "region_high")]
        )[-1]
      )
    }

  }

  if (shade) return(shaded_equivalence_plot(results, ...))

  unshaded_equivalence_plot(results, ...)

}
