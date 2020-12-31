errdist <- function(n, errpar, cn = 0) {
    eps = 0.2
    sigmac = 5
    ind = rbinom(n, 1, eps)
    if (cn == 0) {
        errdist = rnorm(n, errpar[1], errpar[2])
    }
    if (cn == 1) {
        x = rnorm(n, errpar[1], errpar[2])
        errdist = x * (1 - ind) + sigmac * x * ind
    }
    errdist
}
