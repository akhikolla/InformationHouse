#' Classify activity intensity
#'
#' Supports intensity classification via energy expenditure with or without
#' additional posture requirements (i.e., for sedentary behavior to be in
#' lying/seated posture)
#'
#' @param mets numeric vector of metabolic equivalents to classify
#' @param posture character vector of postures
#' @param ... further arguments passed to \code{cut}
#'
#' @details
#' If \code{breaks} and \code{labels} arguments are not provided, default values
#' are <= 1.5 METs for sedentary behavior, 1.51-2.99 METs for light physical
#' activity, and >= 3.0 METs for moderate-to-vigorous physical activity.
#'
#' It is expected for the elements of \code{posture} to be one of \code{c("lie",
#' "sit", "stand", "other")}. The function will run (with a warning) if that
#' requirement is not met, but the output will likely be incorrect.
#'
#' @return a factor giving intensity classifications for each element of
#'   \code{mets}
#' @export
#'
#' @examples
#' mets <- seq(1, 8, 0.2)
#' posture <- rep(
#' c("lie", "sit", "stand", "other"), 9
#' )
#'
#' intensity_no_posture <- get_intensity(mets)
#' intensity_posture <- get_intensity(mets, posture)
#' head(intensity_no_posture)
#' head(intensity_posture)
#'
get_intensity <- function(
  mets, posture = NULL, ...
) {

  dots <- list(...)

  ## Assign `cut` parameter values (using defaults if necessary)

    if ("x" %in% names(dots)) {
      warning("Overwriting value of `x` with `mets`")
    }
    dots$x <- round(mets, 2)

    if (!"breaks" %in% names(dots)) {
      dots$breaks <- c(-Inf, 1.51, 3, Inf)
    }

    if (!"labels" %in% names(dots)) {
      dots$labels <- c("SB", "LPA", "MVPA")
    }

    if (!"right" %in% names(dots)) {
      dots$right <- FALSE
    }

  ## Test `cut` parameter inputs

    stopifnot(
      length(dots$breaks) ==
        length(dots$labels) + 1,
      is.logical(dots$right),
      length(dots$right) == 1
    )

  ## Make initial intensity classification
  ##   (which is final if posture is NULL)

    intensity <- do.call(cut, dots)

    if (is.null(posture)) {
      return(intensity)
    }

  ## Make posture corrections, if applicable

    expected_postures <- c(
      NA, "lie", "sit", "stand", "other"
    )

    posture_test <- posture %in% expected_postures
    if (!all(posture_test)) {
      warning(paste(
        "Expecting all postures labeled as `lie`,",
        "`sit`, `stand`, or `other`.\n  Recode using",
        "`base::switch()` prior to calling `get_intensity`,",
        "\n  or else most (or all) postures will be labeled `other`,",
        "which will\n  preclude SB classification"))
      posture <- ifelse(posture_test, posture, "other")
    }

    stopifnot(
      length(posture) == length(dots$x)
    )

    intensity <- as.character(intensity)
    stopifnot(
      all(c("SB", "LPA") %in% dots$labels)
    )

    sb_postures <- c("sit", "lie")

    for (i in seq(intensity)) {

      if (
        !(is.na(posture[i]) |
          is.na(intensity[i]))
      ) {
        if (
          !posture[i] %in% sb_postures &
            intensity[i] == "SB"
        ) intensity[i] <- "LPA"
      }
    }

    factor(intensity, levels = dots$labels)

}
