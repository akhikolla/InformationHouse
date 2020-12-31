#' Estimate fixed-effect variance for Joint Rank Method (JR) in three-level
#' nested design.
#' 
#' Fixed effect variance estimation for Joint Rank Method (JR). It assumes
#' Compound Symmetric (CS) structure of error terms. For k-level design, there
#' are k-1 intra/inter-class parameters to place in a correlation matrix of
#' errors.
#' 
#' Correlation coefficients are obtained using Moment Estimates. See Klole et.
#' al (2009), Bilgic (2012) and HM (2012) 
#' 
#' @param x Data frame of covariates.
#' @param school A vector of cluster. 
#' @param tauhat This is obtained from Rank-based fitting. 
#' \code{tauhat} here~~
#' @param v1 This is 1, main diagonal element for correlation matrix of
#' observations. Correlation of an observation with itself is 1. 
#' @param v2 Intra-cluster correlation coefficient. 
#' @param v3 Intra-subcluster correlation coefficient. 
#' @param section A vector of subclusters, nx1. 
#' @param mat A matrix of numbers of observations in subclusters.  Dimension is
#' Ixmax(number ofsubclusters). Each row indicates one cluster.  
#' @return \item{var}{ The variance of fixed estimated. }
#' @author Yusuf Bilgic
#' @references Y. K. Bilgic. Rank-based estimation and prediction for mixed
#' effects models in nested designs. 2012. URL
#' http://scholarworks.wmich.edu/dissertations/40. Dissertation.
#' 
#' J. Kloke, J. W. McKean and M. Rashid. Rank-based estimation and associated
#' inferences for linear models with cluster correlated errors. Journal of the
#' American Statistical Association, 104(485):384-390, 2009.
#' 
#' T. P. Hettmansperger and J. W. McKean. Robust Nonparametric Statistical
#' Methods. Chapman Hall, 2012. 
#' 
#' @export
beta_var <- function(x, school, tauhat, v1, v2, v3, section, mat) {
    x <- as.matrix(x)
    ublock <- unique(school)
    I <- length(ublock)
    p <- ncol(x)
    sig2 <- matrix(rep(0, p^2), ncol = p)
    for (i in 1:I) {
        x_i <- as.matrix(x[school == ublock[i], ])
        n_i <- nrow(x_i)
        s_i <- bmat_schC(v1, v2, v3, mat[i,])
        sig2 <- sig2 + t(x_i) %*% s_i %*% x_i
    }
    xx_inv <- solve(crossprod(x))
    V <- tauhat^2 * xx_inv %*% sig2 %*% xx_inv
    list(var = V)
}

Bmat_sch <- function(v1, v2, v3, section) {
    vblock <- unique(section)
    m <- length(vblock)
    N <- length(section)
    B <- matrix(0, nrow = N, ncol = N)
    ind <- 1
    for (k in 1:m) {
        nk <- sum(section == vblock[k])
        if (nk == 1) {
            B[ind:(ind + nk - 1), ind:(ind + nk - 1)] <- 1
        }
        else {
            Bk <- diag(rep(v1 - v2, nk)) + v2
            B[ind:(ind + nk - 1), ind:(ind + nk - 1)] <- Bk
        }
        ind <- ind + nk
    }
    B[B == 0] <- v3
    B
}

interc_se <- function(x, ehat, taus, sec, mat) {
    #rho1 <- rhosect(sign(ehat), school, section)
    
    #c1 <- sum(sum((apply(rho1$rho2, 1, sum, na.rm = T)/apply(rho1$npair, 
    #    1, sum)) * t(apply(rho1$mat, 1, sum))))/sum(sum(rho1$mat))
    
    rho1 <- rhosectC(sign(ehat), sec, mat)
    
    c1 <- sum(sum((apply(rho1$rho2, 1, sum, na.rm = T)/apply(choose(mat, 2), 
      1, sum)) * t(apply(mat, 1, sum))))/sum(sum(mat))
  
    
    #rho2 <- rhosch(sign(ehat), school, section)
    
    rho2 <- rhoschC(sign(ehat), sec, mat)
    
    c2 <- sum(((rho2$rho2/rho2$npair) * rho2$nvec)/sum(rho2$nvec), 
        na.rm = T)
    
    N <- sum(sum(mat))
    a <- N^2 - dim(x)[2]
    var_alpha <- (taus^2/a) * (N + 2 * sum(sum(choose(mat, 2))) * 
        c1 + 2 * sum(rho2$npair) * c2)
    list(var_alpha = var_alpha)
}

jrfit2 <- function(x, y, block) {
    x <- as.matrix(x)
    p <- ncol(x)
    fit <- rfit(y ~ x, symmetric = FALSE)
    betahat <- fit$coef[2:(p + 1)]
    alphahat <- fit$coef[1]
    ehat <- fit$resid
    ahat <- scorewil(ehat)
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
                totals <- total + sign(ek[i]) * sign(ek[j])
            }
        }
    }
    M <- sum(sum(choose(nvec, 2))) - p
    rho <- total/M
    rhos <- totals/M
    taus <- taustar(ehat, p)
    nstar <- sum(sum(nvec * (nvec - 1)))/sum(sum(nvec))
    sigmastar <- 1 + nstar * rhos
    names(alphahat) <- "Intercept"
    coef <- c(alphahat, betahat)
    res <- list(coefficients = coef, residuals = ehat, fitted.values = y - 
        ehat, rhohat = rho, tauhat = fit$tauhat, rhohats = rhos, 
        tauhats = taus, sigmastar = sigmastar, x = x, y = y, 
        block = block)
    class(res) <- list("jrfit")
    res
}

mat_vec <- function(school, section) {
    ublock <- unique(school)
    I <- length(ublock)
    mat <- matrix(0, nrow = I, ncol = max(as.numeric(section)))
    for (i in 1:I) {
        vblock <- unique(section[school == ublock[i]])
        a_i <- school[school == ublock[i]]
        a_i <- a_i[!is.na(a_i)]
        for (j in 1:length(vblock)) {
            a_ij <- a_i[section[school == ublock[i]] == vblock[j]]
            a_ij <- a_ij[!is.na(a_ij)]
            mat[i, j] <- length(a_ij)
        }
    }
    return(mat)
}

pairup1 <- function(x, type = "less") {
    x = as.matrix(x)
    n = dim(x)[1]
    i = rep(1:n, rep(n, n))
    j = rep(1:n, n)
    c1 = apply(x, 2, function(y) {
        rep(y, rep(length(y), length(y)))
    })
    c2 = apply(x, 2, function(y) {
        rep(y, length(y))
    })
    ans = cbind(c1, c2)
    ans = ans[(i < j), ]
    n = dim(ans)[1]
    list(ans = ans, n = n)
}

pairup2 <- function(x, y) {
    x = as.matrix(x)
    y = as.matrix(y)
    n1 = dim(x)[1]
    n2 = dim(y)[1]
    c1 = rep(x, rep(n2, n1))
    c2 = rep(y, n1)
    ans = cbind(as.numeric(c1), as.numeric(c2))
    list(ans = ans, n = n1 * n2)
}


#' Cluster Correlation Coefficient Estimate
#' 
#' Moment estimate version of correlation coefficient in a cluster in a
#' three-level nested design.
#' 
#' 
#' @param ahat A vector of scores. Wilcoxon scores are used in the package.
#' @param school A vector of clusters.
#' @param section A vector of subclusters.
#' @references Y. K. Bilgic. Rank-based estimation and prediction for mixed
#' effects models in nested designs. 2012. URL
#' http://scholarworks.wmich.edu/dissertations/40. Dissertation.
#' 
#' @export 
rhosch <- function(ahat, school, section) {
    ublock <- unique(school)
    I <- length(ublock)
    nvec <- vector(I, mode = "numeric")
    rho1_vec <- rho2_vec <- rho3_vec <- rho4_vec <- npair <- nvec
    for (i in 1:I) {
        vblock <- unique(section[school == ublock[i]])
        
        
        m <- choose(length(vblock), 2)
        nn <- 0
        rho1 <- rho2 <- rho3 <- rho4 <- 0
        a_i <- ahat[school == ublock[i]]
        
        
        for (j in 1:(length(vblock) - 1)) {
            for (k in (j + 1):length(vblock)) {
                if (length(vblock) == 1) {
                  a_ij <- a_i[section[school == ublock[i]] == 
                    vblock[j]]
                  a_ik <- a_ij
                  pair = 1
                  nn <- 0
                  rho1 <- 0
                  rho2 <- 0
                  rho3 <- 0
                  rho4 <- 0
                }
                else {
                  a_ij <- a_i[section[school == ublock[i]] == 
                    vblock[j]]
                  a_ik <- a_i[section[school == ublock[i]] == 
                    vblock[k]]
                  pair = as.matrix(pairup2(a_ij, a_ik)$ans)
                  nn <- nn + pairup2(a_ij, a_ik)$n
                  rho1 <- rho1 + as.numeric(cor(as.numeric(pair[, 
                    1]), as.numeric(pair[, 2])))
                  rho2 <- rho2 + sum(as.numeric(pair[, 1]) * 
                    as.numeric(pair[, 2]))
                  rho3 <- rho3 + (sum(as.numeric(pair[, 1]) * 
                    as.numeric(pair[, 2])) - length(pair[, 1]) * 
                    mean(a_ij) * mean(a_ik))/(sd(pair[, 1] - 
                    mean(a_ij)) * sd(pair[, 2] - mean(a_ik)))
                  rho4 <- rho4 + (sum(pair[, 1] * pair[, 2]) - 
                    length(pair[, 1]) * mean(a_i) * mean(a_i)^2)/(sd(pair[, 
                    1] - mean(a_i)) * sd(pair[, 2] - mean(a_i)))
                }
            }
        }
        rho1_vec[i] = rho1
        rho2_vec[i] = rho2
        rho3_vec[i] = rho3
        rho4_vec[i] = rho4
        nvec[i] <- length(a_i)
        npair[i] = nn
    }
    list(rho1 = rho1_vec, rho2 = rho2_vec, rho3 = rho3_vec, rho4 = rho4_vec, 
        nvec = nvec, npair = npair)
}



#' Subcluster Correlation Coefficient Estimate
#' 
#' Moment estimate version of correlation coefficient in a subcluster in a
#' three-level nested design.
#' 
#' @param ahat A vector of scores. Wilcoxon scores are used in the package.
#' @param school A vector of clusters. 
#' @param section A vector of subclusters.
#' @references Y. K. Bilgic. Rank-based estimation and prediction for mixed
#' effects models in nested designs. 2012. URL
#' http://scholarworks.wmich.edu/dissertations/40. Dissertation.
#' 
#' @export
rhosect <- function(ahat, school, section) {
    ublock <- unique(school)
    I <- length(ublock)
    nvec <- matrix(0, nrow = I, ncol = max(as.numeric(section)))
    rho1_vec <- rho2_vec <- rho3_vec <- rho4_vec <- nvec
    for (i in 1:I) {
        vblock <- unique(section[school == ublock[i]])
        a_i <- ahat[school == ublock[i]]
        a_i <- a_i[!is.na(a_i)]
        for (j in 1:length(vblock)) {
            a_ij <- a_i[section[school == ublock[i]] == vblock[j]]
            a_ij <- a_ij[!is.na(a_ij)]
            nvec[i, j] <- length(a_ij)
            pair = pairup(a_ij)
            if (dim(pair)[2] == 1) {
                rho1_vec[i, j] = 0
                rho2_vec[i, j] = 0
                rho3_vec[i, j] = 0
                rho4_vec[i, j] = 0
            }
            else {
                rho1_vec[i, j] = as.numeric(cor(pair[, 1], pair[, 
                  2]))
                rho2_vec[i, j] = sum(pair[, 1] * pair[, 2])
                rho3_vec[i, j] = (1/length(pair[, 1])) * (sum(pair[, 
                  1] * pair[, 2]) - length(pair[, 1]) * mean(a_ij)^2)/(sd(pair[, 
                  1] - mean(a_ij)) * sd(pair[, 2] - mean(a_ij)))
                rho4_vec[i, j] = (1/length(pair[, 1])) * (sum(pair[, 
                  1] * pair[, 2]) - length(pair[, 1]) * mean(a_i)^2)/(sd(pair[, 
                  1] - mean(a_i)) * sd(pair[, 2] - mean(a_i)))
            }
        }
    }
    list(rho1 = rho1_vec, rho2 = rho2_vec, rho3 = rho3_vec, rho4 = rho4_vec, 
        mat = nvec, npair = choose(nvec, 2))
}

sch_vec <- function(school, section) {
    ublock <- unique(school)
    I <- length(ublock)
    vvec <- as.vector(sec_vec(school, section))
    sctn = rep(0, length(school))
    coll = 0
    for (i in 1:length(vvec)) {
        vblock <- unique(section[school == ublock[i]])
        cols <- rep(0, vvec[i])
        for (j in 1:vvec[i]) {
            zz = section[school == ublock[i]] == vblock[j]
            cols[j] = length(zz[zz == TRUE])
        }
        sctn = rep(1:vvec[i], cols)
        coll = cbind(coll, t(sctn))
    }
    coll = coll[-1]
    return(coll)
}

scorewil <- function(ehat) {
    N = length(ehat)
    u = rank(ehat, ties.method = c("random"))/(N + 1)
    scorewil = sqrt(12) * (u - 0.5)
    list(scorewil = scorewil)
}

sec_sch_vec <- function(y, school1, section1) {
    one = 1 + 0 * y
    schsize = aggregate(one, by = list(y), FUN = sum)[, 2]
    school = factor(rep(1:length(schsize), schsize))
    sss = aggregate(one, by = list(school1, section1), FUN = sum)[, 
        3]
    section = factor(rep(1:length(sss), sss))
    ublock <- unique(school)
    I <- length(ublock)
    vvec <- as.vector(sec_vec(school, section))
    sctn = rep(0, length(school))
    coll = 0
    for (i in 1:length(vvec)) {
        vblock <- unique(section[school == ublock[i]])
        cols <- rep(0, vvec[i])
        for (j in 1:vvec[i]) {
            zz = section[school == ublock[i]] == vblock[j]
            cols[j] = length(zz[zz == TRUE])
        }
        sctn = rep(1:vvec[i], cols)
        coll = cbind(coll, t(sctn))
    }
    coll = coll[-1]
    return(coll)
}

sec_vec <- function(school, section) {
    ublock <- unique(school)
    I <- length(ublock)
    vvec <- matrix(0, nrow = I)
    for (i in 1:I) {
        vblock <- unique(section[school == ublock[i]])
        J <- length(vblock)
        vvec[i] <- J
    }
    return(vvec)
}


#' Wilcoxon estimate for independent linear models
#' 
#' This function gets weighted rank based fittings.
#' 
#' 
#' @param y Response vector of nx1.
#' @param x Design matrix, pxn, without intercept.
#' @references J. T. Terpstra and J. W. McKean. Rank-based analysis of linear
#' models using R. Journal of Statistical Software, 14(7) 1 -- 26, 7 2005. ISSN
#' 1548-7660. URL http://www.jstatsoft.org/v14/i07.
#' 
#' @export
wilonestep <- function(y, x) {
    n = length(y)
    fitw = wwest(x, y, print.tbl = F)
    theta = fitw$tmp1$coef
    ehat = fitw$tmp1$residuals
    taus = taustar(ehat, p = dim(x)[2], conf = 0.95)
    tauhat = wilcoxontau(ehat, p = dim(x)[2])
    list(theta = theta, ehat = ehat, tauhat = tauhat, taus = taus)
}
