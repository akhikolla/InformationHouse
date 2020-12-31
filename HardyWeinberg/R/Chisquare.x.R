Chisquare.x <- function (X) 
{
    mono <- FALSE
    n <- sum(X)
    lab <- names(X)
    nfAA <- X[lab == "AA"]
    nfAB <- X[lab == "AB"]
    nfBB <- X[lab == "BB"]
    nmA <- X[lab == "A"]
    nmB <- X[lab == "B"]
    nm <- nmA + nmB
    nf <- n - nm
    X <- c(nmA, nmB, nfAA, nfAB, nfBB)
    names(X) <- c("A", "B", "AA", "AB", 
                  "BB")
    p <- (2 * nfAA + nfAB + nmA)/(2 * nf + nm)
    theta <- nm/n
    names(p) <- NULL
    q <- 1 - p
    if ((p == 0) | (q == 0)) {
      mono <- TRUE
    }
    obs <- X
    expected <- c(theta * n * p, theta * n * (1 - p), (1 - 
                                                         theta) * n * p^2, (1 - theta) * n * 2 * p * q, (1 - 
                                                                                                           theta) * n * q^2)
    names(expected) <- names(X)
    chi <- (abs(obs - expected))^2
    if (!mono) {
      chi2 <- chi/expected
      chisq <- sum(chi2)
    }
    else {
      chisq <- NA
    }
    return(chisq)
}
