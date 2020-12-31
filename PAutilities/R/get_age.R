#' Calculate age
#'
#' Takes two Date objects and calculates age based on
#' \code{\link[base]{difftime}} (in days) divided by 365.2425 days per year (for
#' age in years) or 30.4375 days per month (for age in months).
#'
#' @param birthdate Date object giving the date of birth
#' @param current_date Date object giving the date from which age is to be
#'   calculated
#' @param units The units in which age should be reported
#'
#' @return Numeric value giving age in the specified units.
#' @export
#'
#' @examples
#' get_age(as.Date("2000-01-01"), Sys.Date(), "years")
get_age <- function(
  birthdate, current_date, units = c("years", "months")
) {

  stopifnot(all(
    c(class(birthdate), class(current_date)) == "Date"
  ))

  units <- match.arg(
    units, c("years", "months", "Error"),
    several.ok = TRUE
  )

  age_days <- as.numeric(
    difftime(current_date, birthdate, units = "days")
  )

  sapply(
    units,
    function(x) {
      switch(
        x,
        "years" = age_days / 365.2425,
        "months" = age_days / 30.4375
      )
    }
  )

}
