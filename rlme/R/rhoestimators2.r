beta_var2 <- function (x, school, tauhat, v1, v2, sec) {
    x <- as.matrix(x)
    ublock <- unique(school)
    I <- length(ublock)
    p <- ncol(x)
    sig2 <- matrix(rep(0, p^2), ncol = p)
    for (i in 1:I) {
        x_i <- as.matrix(x[school == ublock[i], ])
        n_i <- nrow(x_i)
        #s_i <- Bmat2(v1, v2, school[school == ublock[i]])
        s_i <- bmat2C(v1, v2, c(n_i))
        sig2 <- sig2 + s_i[1] * t(x_i) %*% x_i
    }
    xx_inv <- solve(crossprod(x))
    V <- tauhat^2 * xx_inv %*% sig2 %*% xx_inv
    list(var = V)
}

Bmat2 <- function (v1, v2, block) {
    ublock <- unique(block)
    m <- length(ublock)
    N <- length(block)
    B <- matrix(0, nrow = N, ncol = N)
    ind <- 1
    for (k in 1:m) {
        nk <- sum(block == ublock[k])
        if (nk == 1) {
            B[ind:(ind + nk - 1), ind:(ind + nk - 1)] <- 1
        }
        else {
            Bk <- diag(rep(v1 - v2, nk)) + v2
            B[ind:(ind + nk - 1), ind:(ind + nk - 1)] <- Bk
        }
        ind <- ind + nk
    }
    B
}

jrfit2 <- function (pp, ehat, ahat, block) {
    p = pp - 1
    ublock <- unique(block)
    m <- length(ublock)
    nvec <- vector(m, mode = "numeric")
    totals <- total <- 0
    for (k in 1:m) {
        ak <- ahat[block == ublock[k]]
        ek <- ehat[block == ublock[k]]
        nvec[k] <- length(ak)
        for (i in 1:(nvec[k] - 1)) {
            for (j in (i + 1):nvec[k]) {
                total <- total + ak[i] * ak[j]
                totals <- totals + sign(ek[i]) * sign(ek[j])
            }
        }
    }
    M <- sum(choose(nvec, 2)) - p
    rho <- total/M
    rhos <- totals/M
    taus <- taustar(ehat, p)
    nstar <- sum(nvec * (nvec - 1))/sum(nvec)
    sigmastar <- 1 + nstar * rhos
    list(rho = rho, rhos = rhos, taus = taus, nstar = nstar, 
        sigmastar = sigmastar)
}
