fun_ahrf <- function(oy, od, oz, best, tau, alpha, repnum, ...) {
    n <- length(oy)

    oyl <- c(0, oy[1:(n - 1)])
    oyr <- c(oy[2:n], tau)

    bt <- exp(-best)
    bt1 <- bt[1]
    bt2 <- bt[2]

    jh <- 1e-08
    K <- n:1

    b <- as.numeric(best) + cbind(c(jh, 0), -c(jh, 0), c(0, jh), -c(0,
        jh))
    gamma1 <- exp(-matrix(b[1, ], nrow = 1) %x% oz)
    gamma2 <- exp(-matrix(b[2, ], nrow = 1) %x% oz)
    Lambda1 <- apply(od * gamma1/K, 2, cumsum)
    Lambda2 <- apply(od * gamma2/K, 2, cumsum)
    P <- exp(-Lambda2)
    PL <- rbind(1, P[1:(n - 1), ])
    R <- apply(PL * od * gamma1/K, 2, cumsum)/P

    denom <- gamma1 + gamma2 * R
    u1 <- -(oz * od) %*% (gamma1/denom) + oz %*% (R/denom)
    u2 <- -(oz * od) %*% (gamma2 * R/denom) + oz %*% (log(denom/gamma1)/gamma2) -
        oz %*% (R/denom)
    qf <- rbind(u1, u2)
    pq <- cbind(qf[, 1] - qf[, 2], qf[, 3] - qf[, 4])/2/jh
    pr <- cbind(R[, 1] - R[, 2], R[, 3] - R[, 4])/2/jh
    r <- R[, 1]
    rl <- c(0, r[1:(n - 1)])
    dr <- r - rl

    po <- P[, 1]
    plo <- PL[, 1]

    g1 <- gamma1[, 1]
    g2 <- gamma2[, 1]

    den <- bt1 + bt2 * r
    denl <- bt1 + bt2 * rl

    hr <- (1 + r)/den
    ntl <- sum(oy < tau)
    nt <- sum(oy <= tau)

    A <- (1 + r[nt]) * den[nt]
    B <- log(bt1 * hr[nt])/(bt1 - bt2)
    ahrw <- B * (1 + r[nt])/r[nt]
    tem <- bt2 * r[nt]/den[nt]/(bt1 - bt2)
    At1 <- (bt1 * ahrw - bt2 * hr[nt])/(bt1 - bt2)
    At2 <- bt2 * (hr[nt] - ahrw)/(bt1 - bt2)
    C <- (1/den[nt] - B/r[nt])/r[nt]
    B1 <- c(At1, At2) + C * pr[nt, ]
    B2 <- C/po[nt]

    inrs1 <- c()
    inrs2 <- c()
    inrw <- c()
    inrsw1 <- c()
    inrsw2 <- c()

    for (ti in 1:n) {
        yk <- (oy >= oy[ti])
        dk <- g1 + r[ti] * g2
        tek <- yk/dk
        inrs1[ti] <- t(oz) %*% (g1 * tek/dk)
        inrs2[ti] <- R[ti] * t(oz) %*% (g2 * tek/dk)
        inrw[ti] <- t(g2) %*% tek
        inrsw1[ti] <- t(oz) %*% (g1 * g2 * tek/dk^2)
        inrsw2[ti] <- R[ti] * t(oz) %*% (g2^2 * tek/dk^2)
    }

    inr1 <- (inrsw1 - inrw * inrs1/K) * dr/po
    inr2 <- (inrsw2 - inrw * inrs2/K) * dr/po

    inr1 <- inr1 + sum(inr1) - cumsum(inr1)
    inr2 <- inr2 + sum(inr2) - cumsum(inr2)

    rmul <- plo/K
    inr1 <- inr1 * rmul
    inr2 <- inr2 * rmul
    di <- g1 + g2 * r

    inq <- solve(-pq/n)

    xi1d <- oz * g1/di - inrs1 * di/K + inr1 * di
    xi2d <- oz * g2 * r/di - inrs2 * di/K + inr2 * di
    xi1d <- xi1d * od
    xi2d <- xi2d * od

    cid <- plo/K * di * od

    B1inq <-  as.numeric(B1 %*% inq)
    cids1 <- B2 * sqrt(n)*cid[1:nt]
    cids2 <- rep(0, length(cids1))
    mb <- ransamf(repnum = repnum, n = n, B1inq = B1inq, xi1d = xi1d, xi2d = xi2d, cids1 = cids1, cids2 = cids2)

    stmb <- sd(mb)
    ca3 <- qnorm(1 - alpha/2)
    uppc <- exp(ca3 * stmb/sqrt(n)/ahrw) * ahrw
    lowc <- exp(-ca3 * stmb/sqrt(n)/ahrw) * ahrw
    zv <- sqrt(n) * log(ahrw) * ahrw/stmb

    result <- list()
    result$estimate <- ahrw
    result$lower <- lowc
    result$upper <- uppc
    result$z <- as.numeric(zv)
    result$pvalue <- 2 * (1 - pnorm(abs(zv)))
    return(result)
}

