#' Dispatch resting metabolic rate calculations based on data frame input that
#' specifies conversions etc.
#'
#' @param calculation data frame giving instructions and values for the
#'   calculations
#' @param weights list of weights to choose from for the calculation
#'
#' @keywords internal
#'
metabolic_row_wise <- function(calculation, weights) {
  # calculation <- calculations[5, ]
  # weights <- all_weights

  if (calculation$method == "Schofield") {
    select_name <- switch(
      calculation$equation,
      "ht_wt" = "weight_height",
      "wt" = "weight_only"
    )
  } else {
    select_name <- "FAO"
  }

  weights <- weights[select_name]
  stopifnot(length(weights) == 1)

  switch(
    select_name,
    "weight_height" = assemble_wt_ht(calculation, weights),
    "weight_only" = assemble_wt(calculation, weights),
    "FAO" = assemble_fao(calculation, weights)
  )
}

#' @rdname metabolic_row_wise
#' @keywords internal
#'
assemble_wt_ht <- function(calculation, weights) {

  basal_value <-
    (weights$weight_height$weight * calculation$Wt) +
    (weights$weight_height$height * calculation$Ht) +
    (weights$weight_height$intercept)

  basal_VO2_mlkgmin <-
    basal_value * calculation$MJ_conversion / 24 / 60 /
    calculation$kcal_conversion / calculation$Wt * 1000

  data.frame(
    method = calculation$method,
    equation = "Weight_and_Height",
    MJ_conversion_char = calculation$MJ_conversion_char,
    MJ_conversion = calculation$MJ_conversion,
    kcal_table = calculation$kcal_table,
    kcal_conversion = calculation$kcal_conversion,
    basal_value = basal_value,
    basal_units = "MJ/day",
    basal_VO2_mlkgmin = basal_VO2_mlkgmin,
    row.names = NULL,
    stringsAsFactors = FALSE
  )

}

#' @rdname metabolic_row_wise
#' @keywords internal
#'
assemble_wt <- function(calculation, weights) {

    basal_value <-
      (weights$weight_only$weight * calculation$Wt) +
      (weights$weight_only$intercept)

    basal_VO2_mlkgmin <-
      basal_value * calculation$MJ_conversion / 24 / 60 /
      calculation$kcal_conversion / calculation$Wt * 1000

    data.frame(
      method = calculation$method,
      equation = "Weight_Only",
      MJ_conversion_char = calculation$MJ_conversion_char,
      MJ_conversion = calculation$MJ_conversion,
      kcal_table = calculation$kcal_table,
      kcal_conversion = calculation$kcal_conversion,
      basal_value = basal_value,
      basal_units = "MJ/day",
      basal_VO2_mlkgmin = basal_VO2_mlkgmin,
      row.names = NULL,
      stringsAsFactors = FALSE
    )

}

#' @rdname metabolic_row_wise
#' @keywords internal
#'
assemble_fao <- function(calculation, weights) {

  basal_value <-
    (weights$FAO$weight * calculation$Wt) +
    (weights$FAO$intercept)

  basal_VO2_mlkgmin <-
    basal_value / 24 / 60 /
    calculation$kcal_conversion / calculation$Wt * 1000

  data.frame(
    method = calculation$method,
    equation = "FAO",
    MJ_conversion_char = calculation$MJ_conversion_char,
    MJ_conversion = calculation$MJ_conversion,
    kcal_table = calculation$kcal_table,
    kcal_conversion = calculation$kcal_conversion,
    basal_value = basal_value,
    basal_units = "kcal/day",
    basal_VO2_mlkgmin = basal_VO2_mlkgmin,
    row.names = NULL,
    stringsAsFactors = FALSE
  )

}
