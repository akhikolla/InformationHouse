#' Round vector of number to percentages
#'
#' Rounds a vector of numeric values to percentages ensuring that they add up to 100%
#'
#' Returns a vector of numeric values.
#'
#' @param x A numeric vector with non-negative values.
#' @param decimals An integer giving the number of decimals that are used
#' @param ties A string that is either 'random' (the default) or 'last'. Determines how to break ties. Random is random, last prefers to break ties at the last position
#' @return Returns a numeric vector of the same length as x
#' @author Claus Ekstrom \email{claus@@rprimer.dk}
#' @keywords manip
#' @examples
#'
#' f <- c(1,2,1,3,2,1,2,3,1)
#' round_percent(f)
#'
#'
#' @export
round_percent <- function(x, decimals=0L, ties=c("random", "last")) {

    ties <- match.arg(ties)

    ## Do a few sanity checks
    if (!is.numeric(x)) { stop("only works on numeric vectors") }
    if (min(x)<0) { stop("only works on non-negative vectors") }
    if (decimals<0) { stop("number of decimals should be a non-negative integer") }
    decimals <- as.integer(decimals)

    multiplier <- 10^(2+decimals)

    x <- x/sum(x)*multiplier  # Standardize result
    res <- floor(x)           # Find integer bits
    rsum <- sum(res)          # Find out how much we are missing
    if(rsum<multiplier) {
        ## Distribute points based on remainders and a random tie breaker
        tiebreaker <- switch(ties,
                             random = sample(length(x)),
                             last = seq(length(x))
                             )
        o <- order(x%%1, tiebreaker, decreasing=TRUE)
        res[o[1:(multiplier-rsum)]] <- res[o[1:(multiplier-rsum)]]+1
    }
    res/(10^decimals)
}
