#' Obtain preference lists for predicted and actual (reference) activity
#' transitions
#'
#' When predicting activity transitions, the behavior of the predictor is not
#' known a priori. It may predict too many or too few transitions, and its
#' "intent" is also unknown. Therefore, some method is necessary in order to
#' determine which predictions (if any) should be taken to correspond to a
#' reference transition. There should also be a record of false positives and
#' false negatives. The problem is treated as an instance of the college
#' admissions problem, wherein both parties give their preferences for who they
#' would like to be matched with, and a stable arrangement is sought. This
#' function supports the overall goal by assigning the preferences based on the
#' temporal proximity of predicted and actual transitions. Preferences beyond a
#' specified \code{window_size} are not allowed.
#'
#' @inheritParams get_transition_info
#' @inheritParams summary.transition
#'
#' @return A list of matrices giving distance-based preferences for both the
#'   predicted and reference transitions, formatted to pass directly to
#'   \code{\link[matchingMarkets]{hri}}
#' @keywords internal
#'
#' @rdname summary.transition
get_preferences <- function(
  predictions, references, window_size, missing_info
) {

  # Indices of transitions

    ref_i <- which(references == 1)
    pred_i <- which(predictions == 1)

  # Initialize

    prefs <- list(
      window_size = window_size,
      student_reference = references,
      college_prediction = predictions,
      student_reference_i = ref_i,
      college_prediction_i = pred_i,
      student_reference_colnames = ref_i,
      college_prediction_colnames = pred_i,
      false_negative_indices = integer(0),
      false_positive_indices = integer(0)
    )

  prune_prefs(prefs)

}
