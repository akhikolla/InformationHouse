onesampwil <- function (x, test = T, alt = 0, theta0 = 0, alpha = 0.05, maktable = T, 
    plotb = T) {
    n = length(x)
    if (test) {
        as = rank(abs(x - theta0))
        ts = sum(sign(x - theta0) * as)
        zp = (ts - 1)/sqrt((n * (n + 1) * (2 * n + 1))/6)
        zn = (ts + 1)/sqrt((n * (n + 1) * (2 * n + 1))/6)
        if (alt == 1) {
            pval = 1 - pnorm(zp)
            zs = zp
        }
        if (alt == -1) {
            pval = pnorm(zn)
            zs = zn
        }
        if (alt == 0) {
            if (ts >= 0) {
                pval = 2 * (1 - pnorm(zp))
                zs = zp
            }
            else {
                pval = 2 * pnorm(zn)
                zs = zn
            }
        }
    }
    #xpairs = pairupC(x, "leq")
    xpairs = pairup(x, type="leq")
    
    was = (xpairs[, 1] + xpairs[, 2])/2
    est = median(was)
    was = sort(was)
    crit = -qnorm(alpha/2)
    mu = n * (n + 1)/4
    sig = sqrt((n * (n + 1) * (2 * n + 1))/24)
    cut = round(mu - 0.5 - crit * sig)
    if (cut < 0) {
        cut = 0
    }
    lci = was[cut + 1]
    uci = was[(n * (n + 1)/2) - cut]
    tau = sqrt(n/(n - 1)) * ((sqrt(n) * (uci - lci))/(2 * crit))
    if (maktable) {
        if (test) {
            cat("\n")
            cat("Results for the Wilcoxon-Signed-Rank procedure", 
                "\n")
            if (alt == 0) {
                cat("Test of theta =", theta0, " versus theta not equal to ", 
                  theta0, "\n")
            }
            if (alt == 1) {
                cat("Test of theta =", theta0, " versus theta greater than ", 
                  theta0, "\n")
            }
            if (alt == -1) {
                cat("Test of theta =", theta0, " versus theta less than ", 
                  theta0, "\n")
            }
            cat("Test-Stat. is T", ts, "Standardized (z) Test-Stat. is ", 
                zs, "p-vlaue", pval, "\n")
            cat("\n")
        }
        cat("Estimate ", est, " SE is ", tau/sqrt(n), "\n")
        pct = 100 * (1 - alpha)
        cat(pct, "%", "Confidence Interval is ", "(", lci, ",", 
            uci, ")", "\n")
        cat("Estimate of the scale parameter tau", tau, "\n")
    }
    if (plotb) {
        boxplot(x)
        title(main = "Boxplot of Data")
    }
    list(ts = ts, zs = zs, pval = pval, est = est, lci = lci, 
        uci = uci, tau = tau)
}
