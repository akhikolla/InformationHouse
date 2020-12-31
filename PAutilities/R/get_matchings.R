#' Obtain the matchings for predicted and actual activity transitions using the
#' college admissions algorithm
#'
#' @inheritParams prune_prefs
#'
#' @return A data frame giving the relative and absolute indices of the
#'   matchings, based on optimal outcomes for the students (i.e., the actual
#'   transitions)
#' @keywords internal
#'
get_matchings <- function(prefs) {

  if (any(
    !length(prefs$student_reference_prefs),
    !length(prefs$college_prediction_prefs)
  )) {
    result <-
      matrix(nrow = 0, ncol = 7) %>%
      data.frame() %>%
      stats::setNames(c(
        "Prediction", "Reference",
        "Prediction_Index", "Reference_Index",
        "abs_lag", "signed_lag", "rejected"
      ))
    return(result)
  }

  ## Get the matchings
  ## (Use of characters (in newer matchingMarkets versions) causes problems with
  ## indexing, so need to test for presence of characters)

    matchings <-
      matchingMarkets::hri(
        s.prefs = prefs$student_reference_prefs,
        c.prefs = prefs$college_prediction_prefs
      ) %>%
      {.$matchings[
        .$matchings$sOptimal == 1, c("college", "student")
      ]} %>%
      sapply(
        function(x) if(is.character(x)) as.numeric(x) else x,
        simplify = FALSE
      ) %>%
      {do.call(
        data.frame,
        c(., stringsAsFactors = FALSE, row.names = NULL)
      )} %>%
      stats::setNames(c("Prediction", "Reference"))

  matchings$Prediction_Index <-
    matchings$Prediction %>%
    prefs$college_prediction_colnames[.]

  matchings$Reference_Index <-
    matchings$Reference %>%
    prefs$student_reference_colnames[.]

  matchings$abs_lag <- abs(
    matchings$Prediction_Index - matchings$Reference_Index
  )

  matchings$signed_lag <- matchings$Prediction_Index -
    matchings$Reference_Index

  data.frame(
    matchings[order(matchings$Prediction), ],
    row.names = NULL
  ) %>%
  sequence_check(.)

}

#' Label matchings for rejection if the assignment violates sequential nature of
#' the data
#'
#' @param prefs a list assembled in \code{\link{get_transition_info}} prior to
#'   assignment as class \code{transition}
#'
#' @return a list ready for assignment as class \code{transition}
#' @keywords internal
sequence_check <- function(prefs) {

  stopifnot(
    all(diff(prefs$Prediction)>0),
    all(diff(prefs$Prediction_Index)>0)
  )

  matchings <- prefs
  matchings$n_crosses <- n_crossings(matchings)
  matchings$rowname <- seq(nrow(matchings))
  # prefs$matchings$rejected <- FALSE

  while (sum(matchings$n_crosses) != 0) {

    drop_index <- order(
      matchings$n_crosses,
      matchings$abs_lag,
      rev(matchings$rowname),
      decreasing = TRUE
    )[1]

    matchings <- matchings[-drop_index, ]

    matchings$n_crosses <- n_crossings(matchings)

  }

  dropped_rows <- setdiff(
    row.names(prefs), row.names(matchings)
  )

  prefs$rejected <- row.names(prefs) %in% dropped_rows

  return(prefs)

}

#' Determine the number of temporal conflicts among the pairs assigned by the
#' Gale-Shapley algorithm
#'
#' @param matchings data frame with pairing information.
#'
#' @return integer vector specifying the number of conflicts for each pairing
#' @keywords internal
n_crossings <- function(matchings) {

  len <- nrow(matchings)
  crossings <- integer(len)

  for (i in seq(len)) {

    current_index <- matchings$Reference_Index[i]

    prev_indices <- 1:i
    prev_indices <- prev_indices[prev_indices != i]
    if (!length(prev_indices)) {
      prev_crossings <- 0
    } else {
      prev_crossings <- sum(
        current_index <= matchings$Reference_Index[prev_indices]
      )
    }

    subs_indices <- i:len
    subs_indices <- subs_indices[subs_indices != i]
    if (!length(subs_indices)) {
      subs_crossings <- 0
    } else {
      subs_crossings <- sum(
        current_index >= matchings$Reference_Index[subs_indices]
        #>= used intentionally here, even though <= was used previously. The
        #pairing algorithm should work such that the = qualifier is never
        #necessary, but even if it is, I'm willing to double count here, for the
        #purposes of making better decisions about which pairings to reject
      )
    }

    crossings[i] <- prev_crossings + subs_crossings

  }

  crossings

}
