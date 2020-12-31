#' Invoke the Transition Pairing Method
#'
#' @param predictions A dummy-coded vector of predicted transitions (1)
#'   interspersed with non-transitions (0). Logical vectors are coerced to
#'   numeric.
#' @param references A dummy-coded vector of actual (i.e., reference)
#'   transitions (1) interspersed with non-transitions (0). Logical vectors are
#'   coerced to numeric.
#' @param window_size The maximum number of indices that are allowed to separate
#'   a predicted and reference transition, before the two are considered
#'   incompatible
#' @param ... additional arguments passed to or from methods, not currently used
#'
#' @return an object of class \code{transition} that contains necessary
#'   information for evaluating the effectiveness of the predictions.
#' @export
#'
#' @note If the lengths of \code{predictions} and \code{references} differ, a
#'   warning is issued, and the shorter vector will be expanded to match the
#'   length of the longer, using the original relative/proportional positions of
#'   the transitions to determine where they should be placed in the expanded
#'   vector. The relative position could be determined different ways, each
#'   having unique implications for how well aligned \code{predictions} and
#'   \code{references} are. Therefore, while this function is not unusable when
#'   the lengths differ, you should make sure you know what you're doing if you
#'   want to use it that way. The safest solution is to expand the shorter
#'   vector yourself.
#'
#' @seealso \code{\link{summary.transition}}
#'
#' @examples
#' set.seed(100)
#' predictions <- (sample(1:100)%%2)
#' references  <- (sample(1:100)%%2)
#' window_size <- 7
#' get_transition_info(predictions, references, window_size)
get_transition_info <- function(
  predictions, references, window_size = 1, ...
) {

  validate_transition_info_input(predictions, references)

  prefs <- get_preferences(
    predictions, references, window_size
  )

  prefs$matchings <- get_matchings(prefs)

  # Clean up the object

  rejects <- prefs$matchings$rejected

  if (!!length(rejects)) {

    prefs$false_negative_indices <-
      prefs$student_reference_i %>%
      {!. %in% prefs$matchings$Reference_Index[!rejects]} %>%
      prefs$student_reference_i[.]

    prefs$false_positive_indices <-
      prefs$college_prediction_i %>%
      {!. %in% prefs$matchings$Prediction_Index[!rejects]} %>%
      prefs$college_prediction_i[.]

  }

  stopifnot(

    nrow(prefs$matchings[!rejects, ]) +
      length(prefs$false_negative_indices) ==
    length(prefs$student_reference_i),

    nrow(prefs$matchings[!rejects, ]) +
      length(prefs$false_positive_indices) ==
    length(prefs$college_prediction_i)

  )

  names(prefs) <-
    names(prefs) %>%
    sapply(recode_trans_names, USE.NAMES = FALSE)

  structure(prefs, class = "transition")

}
