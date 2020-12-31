#' Compute the area under the curve for two vectors.
#'
#' Compute the area under the curve using linear or natural spline
#' interpolation for two vectors where one corresponds to the x values and the
#' other corresponds to the y values.
#'
#' For linear interpolation the auc function computes the area under the curve
#' using the composite trapezoid rule.  For area under a spline interpolation,
#' auc uses the splinefun function in combination with the integrate to
#' calculate a numerical integral. The auc function can handle unsorted time
#' values, missing observations, ties for the time values, and integrating over
#' part of the area or even outside the area.
#'
#' @param x a numeric vector of x values.
#' @param y a numeric vector of y values of the same length as x.
#' @param from The value from where to start calculating the area under the
#' curve. Defaults to the smallest x value.
#' @param to The value from where to end the calculation of the area under the
#' curve. Defaults to the greatest x value.
#' @param type The type of interpolation. Defaults to "linear" for area under
#' the curve for linear interpolation. The value "spline" results in the area
#' under the natural cubic spline interpolation.
#' @param absolutearea A logical value that determines if negative
#' areas should be added to the total area under the curve.  By
#' default the auc function subtracts areas that have negative y
#' values. Set \code{absolutearea=TRUE} to _add_ the absolute value of the negative areas to the total area.
#' @param subdivisions an integer telling how many subdivisions should be used for integrate (for non-linear approximations)
#' @param \dots additional arguments passed on to approx (for linear approximations). In particular rule can be set to determine how values outside the range of x is handled.
#' @return The value of the area under the curve.
#' @author Claus Ekstrom \email{claus@@rprimer.dk}
#' @seealso \code{\link{approx}}, \code{\link{splinefun}},
#' \code{\link{integrate}}
#' @keywords datagen
#' @examples
#'
#' x <- 1:4
#' y <- c(0, 1, 1, 5)
#' auc(x, y)
#'
#' # AUC from 0 to max(x) where we allow for extrapolation
#' auc(x, y, from=0, rule=2)
#'
#' # Use value 0 to the left
#' auc(x, y, from=0, rule=2, yleft=0)
#'
#' # Use 1/2 to the left
#' auc(x, y, from=0, rule=2, yleft=.5)
#'
#' # Use 1/2 to the left with spline interpolation
#' auc(x, y, from=0, rule=2, yleft=.5)
#'
#'
#' @export auc
auc <-
function(x, y, from = min(x, na.rm=TRUE), to = max(x, na.rm=TRUE), type=c("linear", "spline"), absolutearea=FALSE, subdivisions =100, ...)
{
    type <- match.arg(type)

    # Sanity checks
    stopifnot(length(x) == length(y))
    stopifnot(!is.na(from))

    if (length(unique(x)) < 2)
        return(NA)

    if (type=="linear") {

        ## Default option
        if (absolutearea==FALSE) {
            values <- approx(x, y, xout = sort(unique(c(from, to, x[x > from & x < to]))), ...)
            res <- 0.5 * sum(diff(values$x) * (values$y[-1] + values$y[-length(values$y)]))
        } else { ## Absolute areas
            ## This is done by adding artificial dummy points on the x axis
            o <- order(x)
            ox <- x[o]
            oy <- y[o]

            idx <- which(diff(oy >= 0)!=0)
            newx <- c(x, x[idx] - oy[idx]*(x[idx+1]-x[idx]) / (y[idx+1]-y[idx]))
            newy <- c(y, rep(0, length(idx)))
            values <- approx(newx, newy, xout = sort(unique(c(from, to, newx[newx > from & newx < to]))), ...)
            res <- 0.5 * sum(diff(values$x) * (abs(values$y[-1]) + abs(values$y[-length(values$y)])))
        }
        
    } else { ## If it is not a linear approximation
        if (absolutearea)
            myfunction <- function(z) { abs(splinefun(x, y, method="natural")(z)) }
        else
            myfunction <- splinefun(x, y, method="natural")


        res <- integrate(myfunction, lower=from, upper=to, subdivisions=subdivisions)$value
    }

    res
}


