#' Impute placeholder values into a \code{transition} object
#'
#' If missing values are passed to \code{\link{get_transition_info}}, the
#' prediction and criterion data streams become irregularly spaced. This fixes
#' the spacing issue.
#'
#' @param x Effectively a \code{transition} object (passed before class
#'   assignment)
#'
#' @keywords internal
reconstruct_transitions <- function(x) {

  if (!length(x$missing_cases)) return(x)

  orig_length <- sum(
    length(x$student_reference),
    length(x$missing_cases)
  )

  x$student_reference <-
    seq(orig_length) %in% x$student_reference_colnames %>%
    ifelse(1, 0)

  x$college_prediction <-
    seq(orig_length) %in% x$college_prediction_colnames %>%
    ifelse(1, 0)

  x

}
