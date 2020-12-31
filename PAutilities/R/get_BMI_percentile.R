#' Calculate youth BMI percentile from CDC standards
#'
#' @param weight_kg Weight in kilograms
#' @param height_cm height in centimeters
#' @param age_yrs age in years
#' @param age_mos age in months (optional)
#' @param sex Character scalar indicating participant's sex
#' @param output What should be returned: raw percentile, weight status
#'   classification, or both?
#'
#' @return One of: A numeric scalar giving the BMI percentile (for \code{output
#'   = "percentile"}); a factor scalar giving the weight status (for
#'   \code{output = "classification"}); or a list with the percentile and
#'   classification (for \code{output = "both"}).
#' @export
#'
#' @details If \code{age_mos} is \emph{not} provided, it will be calculated
#'   based on \code{age_yrs}, assuming 365.2425 days per year and 30.4375 days
#'   per month. Depending on how the initial age calculation was made, rounding
#'   error will occur. Thus, use of the \code{\link{get_age}} function is
#'   recommended. If \code{age_mos} \emph{is} provided, \code{age_yrs} can be
#'   passed as \code{NULL}.
#'
#' @references
#' This function was developed with reference to public domain resources
#' provided by the Centers for Disease Control and Prevention. For more
#' information, see:
#'
#' \url{https://www.cdc.gov/obesity/childhood/defining.html}
#'
#' \url{https://www.cdc.gov/healthyweight/assessing/bmi/childrens_bmi/tool_for_schools.html}
#'
#' @examples
#' get_BMI_percentile(39.4, 144.5, 12.35, sex = "M")
get_BMI_percentile <- function(
  weight_kg, height_cm, age_yrs,
  age_mos = NULL, sex = c("M", "F"),
  output = c("percentile", "classification", "both")
) {

  sex <- match.arg(sex, c("M", "F", "Error"), several.ok = FALSE)
  stopifnot(sex %in% c("M", "F"))

  reference <- standards[standards$Sex == sex, ]
  stopifnot(
    !any(duplicated(reference$Age)),
    all(diff(order(reference$Age)) == 1)
  )

  if (is.null(age_mos)) {
    daysold <- age_yrs * 365.2425
    age_mos <- daysold / 30.4375
  }

  BMI <- weight_kg / (height_cm ^ 2) * 10000

  increment <- age_mos - floor(age_mos + 0.5) + 0.5

  greater_index <- max(
    which(reference$Age <= (age_mos + 1))
  )

  lesser_index <- max(
    which(reference$Age <= (age_mos))
  )

  L <- (increment * reference$L[greater_index]) +
    ((1-increment)* reference$L[lesser_index])
  M <- (increment * reference$M[greater_index]) +
    ((1-increment)* reference$M[lesser_index])
  S <- (increment * reference$S[greater_index]) +
    ((1-increment)* reference$S[lesser_index])

  Z_score <- (((BMI/M)^L)-1)/(L*S)

  percentile <- unname(
    floor(stats::pnorm(Z_score) * 1000) / 10
  )
  classification <- cut(
    percentile,
    c(-Inf, 5, 85, 95, Inf),
    c("Underweight", "Normal Weight", "Overweight", "Obese"),
    right = FALSE
  )

  output <- match.arg(output)

  switch(
    output,
    "percentile" = percentile,
    "classification" = classification,
    "both" = stats::setNames(
      list(percentile, classification),
      c("percentile", "classification")
    )
  )

}
