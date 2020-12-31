#' Determine the minimum equivalence zone necessary for establishing equivalence
#' in a paired equivalence test
#'
#' Paired equivalence tests are conducted based on a pre-specified zone of
#' equivalence. Following the test, it is useful to know how minimally small or
#' large the zone would have needed to be in order for the test to indicate
#' equivalence.
#'
#' @param result data frame constructed in
#'   \code{\link{paired_equivalence_wrapper}} that provides information about
#'   the paired equivalence test
#'
#' @keywords internal
equivalent_at <- function(result) {

  if (any(result$region_high != abs(result$region_low))) {
    warning(paste(
      "Asymmetrical equivalence region(s) detected, which violates",
      "code\n  assumptions in `equivalent_at`.",
      "This needs fixing."
    ))
  }

  result <- split(result, seq(nrow(result)))

  result <- lapply(
    result, function(x) {

      x$equivalent_at <- get_absolute_equivalent_at(x)
      if (x$scale == "relative") {
        x$equivalent_at <- x$equivalent_at / x$mean_y
      }

      x

    }
  )

  do.call(rbind, result)

}

#' @rdname equivalent_at
get_absolute_equivalent_at <- function(result) {

  eq_at <- max(
    abs(c(result$CI_low, result$CI_high))
  )

  eq_at + (0.001 * eq_at) ## Equivalence zone needs to
                          ## be slightly larger than `eq_at`
                          ## in order to be equivalent

}
