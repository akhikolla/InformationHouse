secdist <- function(n, secpar, cn = 0) 
{
    eps = 0.2
    sigmac = 5
    ind = rbinom(n, 1, eps)
    if (cn == 0) {
        secdist = rnorm(n, secpar[1], secpar[2])
    }
    if (cn == 1) {
        x = rnorm(n, secpar[1], secpar[2])
        secdist = x * (1 - ind) + sigmac * x * ind
    }
    secdist
}
