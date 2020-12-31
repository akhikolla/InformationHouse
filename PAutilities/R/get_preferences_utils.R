#' Rank preferences for an arbitrary proposer and rejecter, based on distance
#' (i.e., difference) between them
#'
#' @param proposer A vector containing indices of transitions for the proposing
#'   (optimal) party
#' @param rejecter A vector containing indices of transitions for the rejecting
#'   (pessimal) party
#' @inheritParams summary.transition
#'
#' @return A matrix with \code{length(proposer)} columns and
#'   \code{length(rejecter)} rows, giving ordered preferences for the proposer,
#'   based on the absolute difference between the ith proposer index and the
#'   rejecter indices, where differences larger than \code{window_size} are
#'   treated as omitted preferences, i.e., non-possibilities.
#' @keywords internal
#'
get_proposer_rank <- function(proposer, rejecter, window_size) {

  if (any(!length(proposer), !length(rejecter))) return(
    matrix(nrow = 0, ncol = 0)
  )

  ranks <-
    proposer %>%
    sapply(
      function(y) rejecter[order(abs(y - rejecter))],
      simplify = FALSE
    ) %>%
    do.call(cbind, .)

  distances <-
    proposer %>%
    sapply(
      function(y) abs(y - rejecter)[order(abs(y - rejecter))],
      simplify = FALSE
    ) %>%
    do.call(cbind, .)

  ranks[distances > window_size] <- NA

  # Convert to relative ranks

  ranks %>%
  apply(
    2, function(y) ifelse(
      y %in% rejecter, match(y, rejecter), y
    )
  ) %>%
  {if (is.vector(.)) matrix(.,nrow=1) else .}

}

#' @rdname prune_prefs
#' @keywords internal
dim_check <- function(x) {

  if (length(x)==1) {

    return(matrix(
      c(x, rep(NA, 3)),
      nrow = 2, ncol = 2
    ))

  }

  if (nrow(x) == 1) {
    x <- rbind(x, NA)
  }

  if (ncol(x) == 1) {
    x <- cbind(x, NA)
  }

  x

}

#' Account for cases that refuse all matches
#'
#' To run the college admissions algorithm, it is assumed that each
#' student/college has at least one possible match. In the activity transition
#' application, possible matches are restricted to "nearby" cases, i.e. within
#' some specified \code{window_size}. Thus, if there are no matches close by,
#' there are no possible matches, and the case needs to be removed from the
#' analysis in the college admissions algorithm. It needs to be re-inserted
#' afterwards and labeled as a false positive or false negative, depending on
#' whether it was a college (i.e., a predicted transition) or a student (i.e.,
#' an actual transition), respectively.
#'
#' @param prefs an object passed from \code{\link{get_preferences}}
#'
#' @return A pruned list of matrices containing only cases with at least one
#'   possible match
#' @keywords internal
#'
prune_prefs <- function(prefs) {

  # Initialize the preferences

    ref_prefs <- get_proposer_rank(
      prefs$student_reference_i,
      prefs$college_prediction_i,
      prefs$window_size
    )

    pred_prefs <- get_proposer_rank(
      prefs$college_prediction_i,
      prefs$student_reference_i,
      prefs$window_size
    )

  # Test for cases to remove

    ref_test <- apply(
      ref_prefs, 2, function(x) all(is.na(x))
    )

    pred_test <- apply(
      pred_prefs, 2, function(x) all(is.na(x))
    )

  # Perform the pruning

    if (any(ref_test)) {
      prefs$false_negative_indices <-
        prefs$student_reference_i[which(ref_test)]
      if (all(ref_test)) {
        prefs$student_reference_colnames <- integer(0)
      } else {
        prefs$student_reference_colnames <-
          prefs$student_reference_i[-which(ref_test)]
      }
    }

    if (any(pred_test)) {
      prefs$false_positive_indices <-
        prefs$college_prediction_i[which(pred_test)]
      if (all(pred_test)) {
        prefs$college_prediction_colnames <- integer(0)
      } else {
        prefs$college_prediction_colnames <-
          prefs$college_prediction_i[-which(pred_test)]
      }
    }

    if (all(!length(ref_test), !length(pred_test))) {
      prefs$false_negative_indices <- prefs$student_reference_i
      prefs$student_reference_colnames <- integer(0)
      prefs$false_positive_indices <- prefs$college_prediction_i
      prefs$college_prediction_colnames <- integer(0)
    }

  # Get the final preferences

    prefs$student_reference_prefs <-
      get_proposer_rank(
        prefs$student_reference_colnames,
        prefs$college_prediction_colnames,
        prefs$window_size
      ) %>%
      dim_check(.)

    prefs$college_prediction_prefs <-
      get_proposer_rank(
        prefs$college_prediction_colnames,
        prefs$student_reference_colnames,
        prefs$window_size
      ) %>%
      dim_check(.)

  prefs

}
