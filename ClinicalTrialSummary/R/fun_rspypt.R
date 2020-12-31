fun_rspypt <- function(oy, od, oz, best, tau, alpha, repnum, ...) {

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

    inq <- solve(-pq/n)

    r <- R[, 1]
    rl <- c(0, r[1:(n - 1)])
    dr <- r - rl
    po <- P[, 1]
    plo <- PL[, 1]
    g1 <- gamma1[, 1]
    g2 <- gamma2[, 1]
    den <- bt1 + bt2 * r

    lamh2 <- log(1 + bt2/bt1 * r)/bt2
    lamh1 <- log(1 + r)
    fbh1 <- 1/(1 + r)
    fbh2 <- exp(-lamh2)

    dlam1 <- lamh1 - c(0, lamh1[1:(n - 1)])
    dlam2 <- lamh2 - c(0, lamh2[1:(n - 1)])
    df1 <- c(1, fbh1[1:(n - 1)]) - fbh1
    df2 <- c(1, fbh2[1:(n - 1)]) - fbh2

    nk <- sum(oy <= tau)

    rst <- as.numeric(fbh1[1:nk] %*% df2[1:nk])
    rsc <- as.numeric(fbh2[1:nk] %*% df1[1:nk])
    rsp <- as.numeric(rst/rsc)

    sct <- fbh1[1:nk] * fbh2[1:nk]
    dsct <- sct - c(1, sct[1:(nk - 1)])
    dlam12 <- dlam2[1:nk] - rsp * dlam1[1:nk]

    denk <- den[1:nk]
    rk <- r[1:nk]

    Bt1 <- rk/denk
    Bt2 <- lamh2[1:nk] - Bt1
    Bt <- cbind(Bt1, Bt2)
    Bt0 <- 1/denk - rsp/(1 + rk)
    H1t <- 1/denk + 1/(1 + rk)
    H2t <- 1/denk - rsp/(1 + rk)
    H3t <- lamh2[1:nk] - rsp * lamh1[1:nk]

    J1t <- Bt + H1t * pr[1:nk, ]
    J2t <- Bt + H2t * pr[1:nk, ]
    dj2t <- J2t[1:nk, ] - rbind(0, J2t[1:(nk - 1), ])
    dh3t <- H3t[1:nk] - c(0, H3t[1:(nk - 1)])

    tBl <- rbind(0, J2t[1:(nk - 1), ])
    B1 <- -t(dlam12 * sct) %*% J1t + sct[nk] * J2t[nk, ] - t(dsct) %*%
        tBl

    h2p <- H2t/po[1:nk]
    h2pl <- c(exp(best[1]) - rsp, h2p[1:(nk - 1)])
    dh2p <- h2p - h2pl
    B2 <- sct[nk] * h2p[nk]

    inCt1 <- sct/po[1:nk] * H1t * dh3t
    inCt2 <- dsct * h2p
    Ct1 <- sum(inCt1) - cumsum(inCt1) + inCt1
    Ct2 <- sum(inCt2) - cumsum(inCt2) + inCt2
    Ct <- -Ct1 - Ct2

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

    xi1d <- oz * g1/di - inrs1 * di/K + inr1 * di
    xi2d <- oz * g2 * r/di - inrs2 * di/K + inr2 * di

    xi1d <- xi1d * od
    xi2d <- xi2d * od

    cid <- plo/K * di * od

    B1inq <- B1 %*% inq
    cids1 <- B2 * sqrt(n) * cid[1:nk]
    cids2 <- sqrt(n) * Ct * cid[1:nk]
    mb <- ransamf(repnum = repnum, n = n, B1inq = B1inq, xi1d = xi1d, xi2d = xi2d, cids1 = cids1, cids2 = cids2)

    stmb <- sd(mb)/rsc
    ca3 <- qnorm(1 - alpha/2)
    uppc <- exp(ca3 * stmb/sqrt(n)/rsp) * rsp
    lowc <- exp(-ca3 * stmb/sqrt(n)/rsp) * rsp
    zv <- sqrt(n) * log(rsp) * rsp/stmb

    result <- list()
    result$estimate <- rsp
    result$lower <- as.numeric(lowc)
    result$upper <- as.numeric(uppc)
    result$z <- as.numeric(zv)
    result$pvalue <- 2 * (1 - pnorm(abs(zv)))

    return(result)
}

