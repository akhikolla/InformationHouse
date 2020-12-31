schdist <- function(n, schpar, cn = 0) 
{
    eps = 0.2
    sigmac = 5
    ind = rbinom(n, 1, eps)
    if (cn == 0) {
        schdist = rnorm(n, schpar[1], schpar[2])
    }
    if (cn == 1) {
        x = rnorm(n, schpar[1], schpar[2])
        schdist = x * (1 - ind) + sigmac * x * ind
    }
    schdist
}
