#' Compute the age of a person from two dates.
#'
#' Compute the age in years of an individual based on the birth date and another date
#'
#' For linear interpolation the auc function computes the area under the curve
#' using the composite trapezoid rule.  For area under a spline interpolation,
#' auc uses the splinefun function in combination with the integrate to
#' calculate a numerical integral. The auc function can handle unsorted time
#' values, missing observations, ties for the time values, and integrating over
#' part of the area or even outside the area.
#'
#' @param from a vector of dates (birth dates)
#' @param to a vector of current dates
#' @return A vector of ages (in years)
#' @author Claus Ekstrom \email{claus@@rprimer.dk}
#' @seealso \code{\link{as.POSIXlt}}
#' @keywords datagen
#' @examples
#'
#' born <- c("1971-08-18", "2000-02-28", "2001-12-20")
#' check <- c("2016-08-28")
#' age(born, check)
#'
#' @export age
age <- function(from, to) {
    from_lt = as.POSIXlt(from)
    to_lt = as.POSIXlt(to)

    age = to_lt$year - from_lt$year

    if (any(age<0))
        stop("ages cannot be negative")

    ifelse(to_lt$mon < from_lt$mon | (to_lt$mon == from_lt$mon & to_lt$mday < from_lt$mday),  age - 1, age)
}

