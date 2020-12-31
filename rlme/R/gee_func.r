weightf <- function(ehats, ahats, med) {
    w = ahats/(ehats - med)
    w[is.infinite(w)] <- 0
    w[is.nan(w)] <- 0
    w[w == 0] <- max(w)
    list(w = w)
}
