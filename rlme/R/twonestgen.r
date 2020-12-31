twonestgen <- function(I, sec, mat, schdist, schpar, secdist, secpar, errdist, 
    errpar, thetatrue) {
    matre = matrix(c(0, 0), ncol = 2)
    school = c(0)
    class = c(0)
    for (i in 1:I) {
        rsch = schdist(1, schpar)
        for (j in 1:sec[i]) {
            rsec = secdist(1, secpar)
            vecsec = rep(rsec, mat[i, j])
            vecsch = rep(rsch, mat[i, j])
            matt = cbind(vecsch, vecsec)
            matre = rbind(matre, matt)
            for (k in 1:mat[i, j]) {
                school = c(school, i)
                class = c(class, j)
            }
        }
    }
    np1 = length(matre[, 1])
    n = np1 - 1
    matre = matre[2:np1, ]
    school = school[2:np1]
    class = class[2:np1]
    xmat = cbind(rbinom(n, 1, 0.5), rnorm(n))
    errs = errdist(n, errpar)
    ey = xmat %*% thetatrue
    y = ey + matre[, 1] + matre[, 2] + errs
    list(y = y, xmat = xmat, ey = ey, matre = matre, errs = errs, 
        school = school, section = class)
}
