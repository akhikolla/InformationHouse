#' @rdname get_transition_info
#' @keywords internal
validate_transition_info_input <- function(predictions, references) {

  ## Test input classes

    if (is.logical(predictions)) {
      predictions <- as.numeric(predictions)
    }

    if (is.logical(references)) {
      references <- as.numeric(references)
    }

  ## Verify there are transitions

    empty_test <- any(
      all(predictions == 0),
      all(references  == 0)
    )

    if (empty_test) {
      warning(
        "No transitions exist in ",
        "`predictions` and/or `references`",
        call. = FALSE
      )
    }

  ## Check for missing values

    if (anyNA(predictions)) {
      warning(
        "Treating ", sum(is.na(predictions)),
        " missing predictions as non-transitions.",
        call. = FALSE
      )
      predictions[is.na(predictions)] <- 0
    }

    if (anyNA(references)) {
      warning(
        "Treating ", sum(is.na(references)),
        " missing references as non-transitions.",
        call. = FALSE
      )
      references[is.na(references)] <- 0
    }

  ## Check for mismatched lengths

    if (length(predictions) != length(references)) {

      warning(paste(
        "`predictions` and `references` have different lengths.\n",
        "  Expanding the shorter (see note in",
        "?PAutilities::get_transition_info)."
      ), call. = FALSE)

      # impute_trans only affects a vector if its length doesn't match
      # `out_length`

      out_length  <- max(length(predictions), length(references))
      predictions <- impute_trans(predictions, out_length)
      references  <- impute_trans(references, out_length)

    }

  ## Finish up

    stopifnot(
      all(predictions %in% c(0, 1)),
      all(references %in% c(0,1))
    )

    assign("predictions", predictions, parent.frame())
    assign("references", references, parent.frame())

    invisible()

}

#' @rdname get_transition_info
#' @keywords internal
impute_trans <- function(x, out_length) {

  if(length(x) == out_length) return(x)

  (which(x==1)/length(x)) %>%
  {.*out_length} %>%
  round() %>%
  {seq(out_length) %in% .} %>%
  ifelse(1, 0) %T>%
  {stopifnot(sum(.) == sum(x, na.rm = TRUE))} %>%
  return()

}

#' @rdname get_transition_info
#' @keywords internal
recode_trans_names <- function(x) {
  switch(
    x,
    "window_size" = "window_size",
    "student_reference" = "references",
    "college_prediction" = "predictions",
    "student_reference_i" = "reference_transition_indices",
    "college_prediction_i" = "prediction_transition_indices",
    "student_reference_colnames" = "pruned_reference_transition_indices",
    "college_prediction_colnames" = "pruned_prediction_transition_indices",
    "false_negative_indices" = "false_negative_indices",
    "false_positive_indices" = "false_positive_indices",
    "student_reference_prefs" = "reference_preferences",
    "college_prediction_prefs" = "prediction_preferences",
    "matchings" = "matchings"
  )
}
