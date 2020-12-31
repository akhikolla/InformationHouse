#' Determine weight status from body mass index
#'
#' Allows users to determine weight status from body mass index (BMI). The
#' function is designed to classify adult weight status, with default settings
#' yielding weight classes defined by the Centers for Disease Control and
#' Prevention (see reference below). Alternatively, the function can be used as
#' a wrapper for \code{\link{get_BMI_percentile}} to obtain classifications for
#' youth.
#'
#' @usage
#'
#' weight_status(bmi = NULL, breaks = c(-Inf, 18.5, 25, 30, 35, 40, Inf),
#'   labels = c("Underweight", "Normal Weight", "Overweight",
#'   "Obese_1", "Obese_2", "Obese_3"), right = FALSE, youth = FALSE, ...)
#'
#' #get_BMI_percentile(weight_kg, height_cm, age_yrs, age_mos = NULL,
#'   #sex = c("M", "F"), output = c("percentile", "classification", "both"))
#'
#' @param bmi numeric scalar. The participant BMI.
#' @param breaks numeric vector. The boundaries for each weight class; passed to
#'   \code{base::cut}, with warnings if \code{-Inf} and \code{Inf} are not
#'   included in the vector.
#' @param labels character vector. The labels for each weight class; passed to
#'   \code{base::cut}, and should have a length one less than the length of
#'   \code{breaks}
#' @param right logical. See \code{?base::cut}
#' @param youth logical. Use function as a wrapper for
#'   \code{\link{get_BMI_percentile}}?
#' @param ... Arguments passed to \code{\link{get_BMI_percentile}}
#'
#' @return a factor scalar reflecting weight status
#' @export
#'
#' @references
#' \url{https://www.cdc.gov/obesity/adult/defining.html}
#'
#' @examples
#' status <- sapply(17:42, weight_status)
#' head(status)
weight_status <- function(
  bmi = NULL,
  breaks = c(-Inf, 18.5, 25, 30, 35, 40, Inf),
  labels = c(
    "Underweight", "Normal Weight", "Overweight",
    "Obese_1", "Obese_2", "Obese_3"
  ),
  right = FALSE,
  youth = FALSE,
  ...
) {

  if (!-Inf %in% breaks) warning(
    "First element of `breaks` should probably be `-Inf`",
    call. = FALSE
  )

  if (!Inf %in% breaks) warning(
    "Last element of `breaks` should probably be `Inf`",
    call. = FALSE
  )

  if (is.null(bmi)) {

    if (!youth) stop(
      "Cannot classify a null BMI. Did you mean",
      " to set `youth = TRUE`?",
      call. = FALSE
    )

  } else {

    if (bmi < 10) warning(
      "BMI < 10 provided. Recommend checking",
      " calculations/units.",
      call. = FALSE
    )

    if (bmi > 80) warning(
      "BMI > 80 provided. Recommend checking",
      " calculations/units.",
      call. = FALSE
    )

  }

  if (youth) {

    get_BMI_percentile(...)

  } else {

    cut(bmi, breaks, labels, right = FALSE)

  }

}
