#' Retrieve conversion factors from kilocalories to oxygen consumption
#'
#' @inheritParams get_bmr
#'
#' @details RER values are matched to the table entries based on the minimum
#'   absolute difference. If there is a tie, the lower RER is taken.
#'
#' @return numeric vector giving the conversion factor from the specified
#'   table(s)
#' @export
#'
#' @references
#' Peronnet, F., & Massicotte, D. (1991). Table of nonprotein respiratory
#' quotient: an update. \emph{Can J Sport Sci}, 16(1), 23-29.
#'
#' Lusk, G. (1924). Analysis of the oxidation of mixtures of carbohydrate and
#' fat: a correction. \emph{Journal of Biological Chemistry}, 59, 41-42.
#'
#' @examples
#' get_kcal_vo2_conversion(0.85, "both")
#'
get_kcal_vo2_conversion <- function(
  RER, kcal_table = c("Lusk", "Peronnet", "both")
) {

  kcal_table <- match.arg(
    kcal_table, c("Lusk", "Peronnet", "both", "ERROR"), TRUE
  )

  if ("both" %in% kcal_table) kcal_table <- c("Lusk", "Peronnet")

  sapply(
    kcal_table,
    function(x) {
      kcal_table <- switch(
        x,
        "Lusk" = lusk,
        "Peronnet" = peronnet
      )
      index <- which.min(
        abs(kcal_table$Non_Protein_RQ  - RER)
      )
      kcal_table$kCal_per_Liter_O2_Uptake[index]
    }
  )

}
