#' @rdname paired_equivalence_test
#' @export
paired_equivalence_test.default <- function(
  x, y, y_type = c("both", "criterion", "comparison"),
  alpha = 0.05, na.rm = TRUE,
  scale = c("relative", "absolute"),
  absolute_region_width = NULL,
  relative_region_width = NULL,
  ...
) {

  stopifnot(is.numeric(x))
  keep <- stats::complete.cases(x, y)

  x <- x[keep]
  y <- y[keep]

  object <- new_paired_equivalence(
    x, y, y_type, alpha, na.rm, scale,
    absolute_region_width, relative_region_width
  )

  if (
    "criterion" %in% object$y_type &&
      object$scale == "absolute"
  ) {

    criterion_results <- paired_equivalence_wrapper(
      object$x, object$y, "criterion",
      object$absolute_region_width,
      object$absolute_region_width,
      object$alpha, object$scale, object$na.rm
    )

  }

  if (
    "criterion" %in% object$y_type &&
      object$scale == "relative"
  ) {

    criterion_results <- paired_equivalence_wrapper(
      object$x, object$y, "criterion",
      object$relative_region_width,
      object$relative_region_width * mean(object$y, na.rm = na.rm),
      object$alpha, object$scale, object$na.rm
    )

  }

  if (
    "comparison" %in% object$y_type &&
      object$scale == "absolute"
  ) {

    comparison_results <- paired_equivalence_wrapper(
      object$x, object$y, "comparison",
      object$absolute_region_width,
      object$absolute_region_width,
      object$alpha, object$scale, object$na.rm
    )

  }

  if ("comparison" %in% y_type && scale == "relative") {

    comparison_results <- paired_equivalence_wrapper(
      object$x, object$y, "comparison",
      object$relative_region_width,
      object$relative_region_width * mean(object$y, na.rm = na.rm),
      object$alpha, object$scale, object$na.rm, FALSE
    )

    ## Do the lower bound test

    lower <- 1 - relative_region_width

    lower_diffs <- x - (lower * y)
    lower_result <- stats::t.test(
      lower_diffs, alternative = "greater", mu = 0
    )$p.value

    ## Do the upper bound test

    upper <- 1/lower

    upper_diffs <- x - (upper * y)
    upper_result <- stats::t.test(
      upper_diffs, alternative = "less", mu = 0
    )$p.value

    comparison_results$tost_p <- max(
      lower_result, upper_result
    )

    comparison_results$tost_sig <- ifelse(
      comparison_results$tost_p < alpha,
      "*", "NS"
    )

  }

  ## Finish up

  results <- rbind(
    criterion_results, comparison_results
  )

  object$results <- results

  object

}
