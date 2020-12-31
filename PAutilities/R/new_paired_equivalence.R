#' A combined constructor and validator for \code{paired_equivalence} objects
#'
#' @inheritParams paired_equivalence_test
#'
#' @return an initialized \code{paired_equivalence} object
#' @keywords internal
new_paired_equivalence <- function(
  x, y, y_type = c("both", "criterion", "comparison"),
  alpha = 0.05, na.rm = TRUE,
  scale = c("relative", "absolute"),
  absolute_region_width = NULL,
  relative_region_width = NULL
) {

  ## Setup for y_type and scale

    y_type <- match.arg(
      y_type,
      c("both", "criterion", "comparison", "Error"),
      TRUE
    )

    if ("both" %in% y_type) {
      y_type <- c("criterion", "comparison")
    }

    scale <- match.arg(scale)

  ## Setup for absolute_region_width

    if (is.null(absolute_region_width) && scale == "absolute") {
      warning(paste(
        "\nAssigning default value (5) for unspecified",
        "width of equivalence region.\n  See `?paired_equivalence_test`."
      ))
      absolute_region_width <- 5
    }

    ## Setup for relative_region_width

    if (is.null(relative_region_width) && scale == "relative") {
      warning(paste(
        "\nAssigning default value (0.1) for unspecified",
        "width of equivalence region.\n  See `?paired_equivalence_test`."
      ))
      relative_region_width <- 0.1
    }

  ## Check for numeric input

    stopifnot(
      is.numeric(y),
      (
        is.numeric(absolute_region_width) |
          is.null(absolute_region_width)
      ),
      (
        is.numeric(relative_region_width) |
          is.null(relative_region_width)
      )
    )

  ## Check for potentially bad specification of relative region

    if (is.numeric(relative_region_width)) {

      stopifnot(relative_region_width >= 0)

      if (relative_region_width > 1) warning(paste(
        "Expecting a relative region width between 0 and 1.",
        "\n  Have you passed a percentage instead of a proportion?"
      ))

    }

  ## Finish up

    object <- list(
      x = x, y = y, alpha = alpha, y_type = y_type,
      na.rm = na.rm,
      absolute_region_width = absolute_region_width,
      relative_region_width = relative_region_width,
      scale = scale
    )

    structure(
      object, class = "paired_equivalence"
    )

}
