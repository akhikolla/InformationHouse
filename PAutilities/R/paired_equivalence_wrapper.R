#' A template function for conducting a paired equivalence test
#'
#' @inheritParams paired_equivalence_test
#' @param width the user-specified width of the equivalence region, possibly a
#'   proportion of an initially-uncalculated mean
#' @param epsilon the calculated width of the equivalence region
#' @param do_test logical. Complete the test? Enables premature return that is
#'   useful in some cases
#'
#' @keywords internal
paired_equivalence_wrapper <- function(
  x, y, y_type, width, epsilon, alpha, scale, na.rm,
  do_test = TRUE
) {

  result <- stats::setNames(

    data.frame(
      matrix(ncol = 14),
      row.names = NULL,
      stringsAsFactors = FALSE
    ),

    c(
      "mean_x", "mean_y", "y_type", "mean_diff", "scale",
      "region_width", "region_low", "region_high", "CI_low",
      "CI_high", "tost_p", "tost_sig", "CI_sig", "equivalent_at"
    )

  )

  result$mean_x <- mean(x, na.rm = na.rm)
  result$mean_y <- mean(y, na.rm = na.rm)
  result$y_type <- y_type
  result$mean_diff <- mean(x - y, na.rm = na.rm)
  result$scale <- scale
  result$region_width <- width

  region_low <- -abs(epsilon)
  region_high <- abs(epsilon)

  result$region_low <- region_low
  result$region_high <- region_high

  if (!do_test) return(result)

  tost_result <- equivalence::tost(
    x = x, y = y, epsilon = epsilon,
    paired = TRUE, var.equal = TRUE,
    conf.level = 1 - alpha, alpha = alpha
  )

  CI_sig <- ifelse(
    (
      tost_result$tost.interval[1] >= region_low &
        tost_result$tost.interval[2] <= region_high
    ),
    "*", "NS"
  )

  result$CI_low <- tost_result$tost.interval[1]
  result$CI_high <- tost_result$tost.interval[2]
  result$tost_p <- tost_result$tost.p.value
  result$tost_sig <- ifelse(
    tost_result$tost.p.value < alpha, "*", "NS"
  )
  result$CI_sig <- CI_sig

  result <- equivalent_at(result)

  result

}
