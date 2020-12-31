combineC <- function (Xmat, alpha = 0.05, fisher = FALSE, varest = 1) 
{
    fvec <- NULL
    varvec <- NULL
    m <- nrow(Xmat)
    for (i in 1:nrow(Xmat)) {
        fhat <- HWf(Xmat[i, ])
        pa <- af(Xmat[i, ])
        n <- sum(Xmat[i, ])
        fvec <- c(fvec, fhat)
        if(varest == 1) varf <- ((1-fhat)^2)*(1-2*fhat)/n + fhat*(1-fhat)*(2-fhat)/(2*n*pa*(1-pa)) else
        stop("unknown option for parameter varest")
        varvec <- c(varvec, varf)
    }
    if (fisher) {
        fvec <- fisherz(fvec)
        varvec <- rep(1/(n - 3), length(fvec))
    }
    MeanTheta <- mean(fvec)
    W <- mean(varvec)
    B <- var(fvec)
    T <- W + (m + 1) * B/m
    gamma <- (1 + 1/m) * B/T
    v <- (m - 1) * (1 + m * W/((m + 1) * B))^2
    stat <- abs(MeanTheta)/sqrt(T)
    pvalimp <- 2 * pt(stat, v, lower.tail = FALSE)
    r <- (1 + 1/m) * B/W
    lambda <- (r + 2/(v + 3))/(r + 1)
    gamma <- (1 + 1/m) * B/T
    gammaalt <- (r + 2/(v + 3))/(r + 1)
    fhatimp <- MeanTheta
    if (fisher) 
        fhatimp <- ifisherz(MeanTheta)
    llf <- fhatimp - qt(1 - alpha/2, v) * sqrt(T)
    ulf <- fhatimp + qt(1 - alpha/2, v) * sqrt(T)
    if (fisher) {
        llf <- fhatimp - qnorm(1 - alpha/2) * sqrt(1/(n - 3))
        ulf <- fhatimp + qnorm(1 - alpha/2) * sqrt(1/(n - 3))
        llf <- ifisherz(llf)
        ulf <- ifisherz(ulf)
    }
    return(list(fhatimp = fhatimp, pvalimp = pvalimp, r = r, 
        lambda = lambda, llf = llf, ulf = ulf, fvec = fvec, varvec = varvec, 
        gammaalt = gammaalt))
}


