#' GEER: General Estimating Equation Rank-Based Estimation Method
#' 
#' The package rlme calls this function for gee method, one of the methods
#' proposed in Bilgic's study (2012). Also see Kloke et al. (2013).  %% ~~ A
#' concise (1-5 lines) description of what the function does. ~~
#' 
#' 
#' @param x Design matrix, pxn, without intercept.
#' @param y Response vector of nx1.
#' @param I Number of clusters.
#' @param sec A vector of subcluster numbers in clusters.
#' @param mat A matrix of numbers of observations in subclusters.  Dimension is
#' Ixmax(number ofsubclusters). Each row indicates one cluster.
#' @param school A vector of clusters, nx1.
#' @param section A vector of subclusters, nx1.
#' @param weight When weight="hbr", it uses hbr weights in GEE weights. By
#'   default, ="wil", it uses Wilcoxon weights. See the theory in the references.
#' @param rprpair By default, it uses "hl-disp" in the random prediction
#'   procedure (RPP). Also, "med-mad" would be an alternative.
#' @param verbose Boolean indicating whether to print out diagnostic messages.
#' 
#' @return
#'   \item{theta}{ Fixed effect estimates. }
#'   \item{ses}{ Standard error for the fixed esimates. }
#'   \item{sigma}{ Variances of cluster, subcluster, and residual. }
#'   \item{ehat}{ Raw error. }
#'   \item{ehats}{ Independence error from last weighted step. }
#'   \item{effect_sch}{ Cluster random error. }
#'   \item{effect_sec}{ Subcluster random error. }
#'   \item{effect_err}{ Epsilon error. }
#'   
#' @author Yusuf K. Bilgic, yekabe@hotmail.com
#' 
#' @seealso rlme, GR_est, JR_est, rprmeddisp
#' 
#' @references Y. K. Bilgic. Rank-based estimation and prediction for mixed
#' effects models in nested designs. 2012. URL
#' http://scholarworks.wmich.edu/dissertations/40. Dissertation.
#' 
#' A. Abebe, J. W. McKean, J. D. Kloke and Y. K. Bilgic. Iterated reweighted
#' rank-based estimates for gee models. 2013. Submitted.
#' 
#' @keywords models
#' 
#' @importFrom magic adiag
#' 
#' @examples
#' 
#' # See the rlme function.
#' 
#' @export 
GEER_est <- function(x, y, I, sec, mat, school, section, weight = "wil", rprpair = "hl-disp", verbose=FALSE) {
    location = 2
    tol <- 1e-23
    weight = tolower(weight)
    if (weight == "wil") {
        weight = 1
    }
    if (weight == "hbr") {
        weight = 2
    }
    
    n = sum(mat)

    fit = minimize_dispersion(as.matrix(x), as.matrix(y))
    b0 = fit$theta
    
    ehat0 = y - cbind(1, x) %*% b0
    fitvc = rprmeddis(I, sec, mat, ehat = ehat0, location, scale, rprpair = rprpair)
    sigmaa2 = fitvc$siga2
    sigmaw2 = fitvc$sigw2
    sigmae2 = fitvc$sigmae2
    x = cbind(1, x)
    collb <- b0
    collsigma <- c(sigmaa2, sigmaw2, sigmae2)
    iter <- 0
    chk <- 1
    b2 <- b0
    max.iter <- 2
    ww <- 0
    
    if(verbose == TRUE) {
      cat("GEER: Starting iterative step\n")
    }
    
    while ((chk > 1e-04) && (iter < max.iter)) {
        iter <- iter + 1
        b1 <- b2
        sigmay = sigymake(I, sec, mat, sigmaa2, sigmaw2, sigmae2)
        sigma12inv = sigmay$sigy12i
        siggma12 <- sigmay$siggma
        
        ystar = sigma12inv %*% y
        xstar = sigma12inv %*% x
        ehats = ystar - xstar %*% b2
        med <- median(ehats)
        ahats = scorewil(ehats)$scorewil
        yss = y - siggma12 %*% rep(1, n) * med
        ys = ystar - siggma12 %*% rep(1, n) * med
        if (weight == 1) {
            w <- weightf(ehats, ahats, med)$w
        }
        if (weight == 2) {
            w <- hbrwts_gr(xstar, ehats)
        }
        ww <- rbind(ww, t(w))
        w <- diag(as.vector(w))
        r <- sigma12inv %*% w %*% sigma12inv
        b2 <- solve(t(x) %*% r %*% x, tol = tol) %*% (t(x) %*% 
            r %*% yss)
        ehat = y - x %*% b2
        fitvc <- rprmeddis(I, sec, mat, ehat = ehat, location, 
            scale, rprpair = rprpair)
        sigmaa2 = fitvc$siga2
        sigmaw2 = fitvc$sigw2
        sigmae2 = fitvc$sigmae2
        chk <- sum((b2 - b1)^2)/sum(b1^2)
        collb <- rbind(collb, c(b2))
        collsigma <- rbind(collsigma, c(sigmaa2, sigmaw2, sigmae2))
    }
    iter = iter + 1
    b <- collb[iter, ]
    theta = b
    ehat = y - x %*% b
    ahat = scorewil(ehat)$scorewil
    tauhat = wilcoxontau(ehat, p = dim(x)[2], verbose=verbose)
    fitvc <- rprmeddis(I, sec, mat, ehat = ehat, location, scale, 
        rprpair = rprpair)
    sigmaa2 = fitvc$siga2
    sigmaw2 = fitvc$sigw2
    sigmae2 = fitvc$sigmae2
    sigma <- c(sigmaa2, sigmaw2, sigmae2)
    effect_sch = fitvc$frei
    effect_sec = fitvc$frew
    effect_err = fitvc$free
    
    sigmay = sigymake(I, sec, mat, sigmaa2, sigmaw2, sigmae2)
    
    sigma12inv = sigmay$sigy12i
    siggma12 = sigmay$siggma
    
    
    if (weight == 1) {
        w <- weightf(ehats, ahats, med)$w
    }
    if (weight == 2) {
        w <- hbrwts_gr(xstar, ehats)
    }
    w <- diag(as.vector(w))
    r <- sigma12inv %*% w %*% sigma12inv
    ystar = sigma12inv %*% y
    xstar = sigma12inv %*% x
    ehats = ystar - xstar %*% b
    ahats = scorewil(ehats)$scorewil
    tauhats = wilcoxontau(ehats, p = dim(x)[2])
    m01 = solve(t(x) %*% r %*% x, tol = tol)
    m03 = solve(t(x) %*% sigma12inv %*% sigma12inv %*% x, tol = tol)
    sisi3 = diag(rep(1, n))
    m15 = (t(x) %*% sigma12inv %*% sisi3 %*% sigma12inv %*% x)
    varb = tauhats^2 * m03 %*% m15 %*% m03
    se31 = sqrt(diag(varb))
    rho1 = rhosect(ahats, school, section)
    rho1_est_9 = sum(sum((apply(rho1$rho2, 1, sum)/apply(rho1$npair, 
        1, sum)) * t(apply(rho1$mat, 1, sum))))/sum(sum(rho1$mat))
    
    #rho2 = rhosch(ahats, school, section)
    rho2 = rhoschC(ahats, sec, mat)
    
    aaa = ((rho2$rho2/rho2$npair) * rho2$nvec)
    aaa = aaa[!is.na(aaa)]
    rho2_est_5 = sum(aaa/sum(rho2$nvec))
    v1 = 1
    v2 = rho1_est_9
    v3 = rho2_est_5
    sisi4 = 0
    
    for (i in unique(school)) {
        sisi4 = adiag(sisi4, bmat_schC(v1, v2, v3, mat[i,]))
    }
    sisi4 = sisi4[2:dim(sisi4)[1], 2:dim(sisi4)[2]]
    sisi5 = matrix(sisi4, ncol = n)
    m16 = (t(x) %*% sigma12inv %*% sisi4 %*% sigma12inv %*% x)
    
    if (weight == 2) {
        varb = m01 %*% m16 %*% m01
        se41 = sqrt(diag(varb))
    }
    if (weight == 1) {
        varb = tauhats^2 * m03 %*% m16 %*% m03
        se41 = sqrt(diag(varb))
    }
    
    list(theta = theta,
         ses_AP = se31,
         ses_CS = se41,
         varb = varb, 
         sigma = sigma,
         ehat = ehat,
         effect_sch = effect_sch, 
         effect_sec = effect_sec,
         effect_err = effect_err,
         iter = iter, 
         w = diag(w))
}



#' GR Method
#' 
#' Fits a model using the GR method
#' 
#' 
#' @param x Covariate matrix or data frame.
#' @param y Response matrix or data frame.
#' @param I Number of clusters 
#' @param sec A vector of subcluster numbers in clusters.
#' @param mat A matrix of numbers of observations in subclusters.  Dimension is
#' Ixmax(number ofsubclusters). Each row indicates one cluster.
#' @param school A vector of clusters, nx1.
#' @param section A vector of subclusters, nx1.
#' @param rprpair By default, it uses "hl-disp" in the random prediction
#' procedure (RPP). Also, "med-mad" would be an alternative.
#' @param verbose Boolean indicating whether to print out messages from the
#' algorithm.
#' 
#' @return
#'   \item{theta}{ Fixed effect estimates. }
#'   \item{ses}{ Standard error for the fixed esimates. }
#'   \item{sigma}{ Variances of cluster, subcluster, and residual. }
#'   \item{ehat}{ Raw error. }
#'   \item{ehats}{ Independence error from last weighted step. }
#'   \item{effect_sch}{ Cluster random error. }
#'   \item{effect_sec}{ Subcluster random error. }
#'   \item{effect_err}{ Epsilon error. }
#'   
#' @author Yusuf Bilgic
#' @keywords models
#' @examples
#' 
#' # See rlme function
#' 
#' @export
GR_est <- function(x, y, I, sec, mat, school, section, rprpair = "hl-disp", verbose=FALSE) {
  
    init = T
    sigmaa2 = 1
    sigmaw2 = 1
    sigmae2 = 1
    thetaold = c(0)
    
    if(nrow(x) < 3000) {
      numstp = 10
    } else {
      numstp = 2
    }
    eps = 1e-04
    iflag2 = 0
    i = 0
    is = 0
    J = sum(sec)
    n = length(y)
    if (is.null(dim(x)[2])) {
        nopar_scale = 1 + 1 + 3
    }
    else {
        nopar_scale = dim(x)[2] + 1 + 3
    }
    coll = matrix(rep(0, nopar_scale), ncol = nopar_scale)
    collresch = matrix(rep(0, I), ncol = I)
    collresec = matrix(rep(0, J), ncol = J)
    collreserr = matrix(rep(0, n), ncol = n)
    while (is == 0) {
        i = i + 1
        
        if(verbose == TRUE) {
          cat(paste0("GR: iteration ", i, "\n"))
        }
        
        if (i == 1) {
            if(verbose == TRUE) {
              cat("GR: getting initial fit from wilstep\n")
            }
            
            wilfit = wilstep(I, sec, mat, init, y, x, sigmaa2, 
                sigmaw2, sigmae2, thetaold, eps, iflag2, rprpair = rprpair)
            coll = rbind(coll, c(wilfit$theta, wilfit$sigmaa2, 
                wilfit$sigmaw2, wilfit$sigmae2))
            effect_sch = wilfit$rea
            effect_sec = wilfit$rew
            effect_err = wilfit$ree
            thetaold = wilfit$theta
            ehat = wilfit$ehat
            theta = coll[dim(coll)[1], 1:(dim(x)[2] + 1)]
            sigma = c(wilfit$sigmaa2, wilfit$sigmaw2, wilfit$sigmae2)
        }
        
        sigmaa2 = wilfit$sigmaa2
        sigmaw2 = wilfit$sigmaw2
        sigmae2 = wilfit$sigmae2
        
        init = F
        
        if(verbose == TRUE) {
          cat("GR: getting step from wilstep\n")
        }
        
        wilfit = wilstep(I, sec, mat, init, y, x, sigmaa2, sigmaw2, 
            sigmae2, thetaold, eps, iflag2, rprpair = rprpair)
        
        if ((wilfit$iflag2 == 1) | (i > numstp)) {
            is = 1
        }
        coll = rbind(coll, c(wilfit$theta, wilfit$sigmaa2, wilfit$sigmaw2, 
            wilfit$sigmae2))
        effect_sch = wilfit$rea
        effect_sec = wilfit$rew
        effect_err = wilfit$ree
        thetaold = wilfit$theta
        ehat = wilfit$ehat
        theta = coll[dim(coll)[1], 1:(dim(x)[2] + 1)]
        sigma = c(wilfit$sigmaa2, wilfit$sigmaw2, wilfit$sigmae2)
        i
    }
    
    if(verbose == TRUE) {
      cat("GR: done iterating\n");
    }
    
    
    if(verbose == TRUE) {
      cat("GR: building sigma\n")
    }
    
    sigmay = sigymake(I, sec, mat, sigma[1], sigma[2], sigma[3])
    sigma12inv = sigmay$sigy12i
    xstar = sigma12inv %*% cbind(1, x)
    ystar = sigma12inv %*% y
    ehats = ystar - xstar %*% theta
    ahats = scorewil(ehats)$scorewil
    ehat = y - cbind(1, x) %*% theta
    
    taus = taustar(ehats, p = dim(x)[2], conf = 0.95)
    tauhat = wilcoxontau(ehats, p = dim(x)[2], verbose=verbose)
    
    XXinv <- solve(crossprod(xstar))
    
    x_c = xstar - apply(xstar, 2, mean)
    
    aa = taus^2 * (XXinv %*% t(xstar)) %*% {
        rep(1, n) %*% solve(t(rep(1, n)) %*% rep(1, n)) %*% rep(1, 
            n)
    } %*% (xstar %*% XXinv)
    
    bb = tauhat^2 * XXinv %*% t(xstar) %*% {
        x_c %*% solve(t(x_c) %*% x_c) %*% t(x_c)
    } %*% (xstar %*% XXinv)
    
    varb = aa + bb
    
    ses <- sqrt(diag(varb))
    
    list(theta = theta, 
         ses = ses,
         sigma = sigma,
         varb = varb, 
         ehat = ehat,
         ehats = ehats,
         effect_sch = effect_sch, 
         effect_sec = effect_sec,
         effect_err = effect_err,
         iter = i, 
         coll = coll,
         xstar = xstar,
         ystar = ystar)
}



#' JR Method
#' 
#' Fit a model using the JR method
#' 
#' 
#' @param x Covariate matrix or data frame
#' @param y Response matrix or data frame
#' @param I Number of clusters.
#' @param sec A vector of subcluster numbers in clusters.
#' @param mat A matrix of numbers of observations in subclusters.  Dimension is
#' Ixmax(number ofsubclusters). Each row indicates one cluster.
#' \code{mat} here~~
#' @param school A vector of clusters, nx1.
#' @param section A vector of subclusters, nx1.
#' @param rprpair By default, it uses "hl-disp" in the random prediction
#' procedure (RPP). Also, "med-mad" would be an alternative. 
#' @param verbose Boolean indicating whether to print out diagnostic messages.
#' @return
#'   \item{theta}{ Fixed effect estimates. }
#'   \item{ses}{ Standard error for the fixed esimates. }
#'   \item{sigma}{ Covariate variance estimates using RPP (Groggel and Dubnicka's procedure). }
#'   \item{ehat}{ Raw error. }
#'   \item{effect_sch}{ Cluster random error. }
#'   \item{effect_sec}{ Subcluster random error. }
#'   \item{effect_err}{ Epsilon error. }
#'   
#' @author Yusuf Bilgic
#' @seealso rlme
#' @export
JR_est <- function(x, y, I, sec, mat, school, section, rprpair = "hl-disp", verbose=FALSE) {
    location = 2
    tol <- 1e-23
    pp <- dim(x)[2] + 1
    
    # Fit the fixed effect coefficients
    
    if(verbose == TRUE) {
      cat("JR: minimizing dispersion function")
    }
    
    wilfit = minimize_dispersion(as.matrix(x), as.matrix(y), method='L-BFGS-B', verbose=verbose)
    
    theta = wilfit$theta[1:(pp)]
    ehat = wilfit$ehat
    ahat = scorewil(ehat)$scorewil
    taus = wilfit$taus
    tauhat = wilfit$tauhat
    
    if(verbose == TRUE) {
      cat("JR: calling RPP algorithm")
    }
    
    scale.fit = rprmeddis(I, sec, mat, ehat, location, scale, 
        rprpair = rprpair)
    sigma = c(scale.fit$siga2, scale.fit$sigw2, scale.fit$sigmae2)
    effect_sch = scale.fit$frei
    effect_sec = scale.fit$frew
    effect_err = scale.fit$free
    var_alpha = interc_se(x, ehat, taus, sec, mat)$var_alpha
    
    #rho1 <- rhosect(ahat, school, section)
    #rho1_est_9 <- sum(sum((apply(rho1$rho2, 1, sum, na.rm = T)/apply(rho1$npair, 
    #    1, sum)) * t(apply(rho1$mat, 1, sum))))/sum(sum(rho1$mat))
    
    if(verbose == TRUE) {
      cat("JR: estimating correlations\n")
    }
    
    if(verbose == TRUE) {
      cat("JR: calling rhosect\n")
    }
    
    rho1 <- rhosectC(ahat, sec, mat)
    rho1_est_9 <- sum(sum((apply(rho1$rho2, 1, sum, na.rm = T)/apply(choose(mat, 2), 1, sum)) * t(apply(mat, 1, sum))))/sum(sum(mat))
    
    
    if(verbose == TRUE) {
      cat("JR: calling rhosch\n")
    }
    
    #rho2 <- rhosch(ahat, school, section)
    rho2 <- rhoschC(ahat, sec, mat)
    
    rho2_est_5 <- sum(((rho2$rho2/rho2$npair) * rho2$nvec)/sum(rho2$nvec), na.rm = T)
    v1 = 1
    v2 = rho1_est_9
    v3 = rho2_est_5
    
    if(verbose == TRUE) {
      cat("JR: calling beta_var\n")
    }
    
    V = beta_var(x, school, tauhat, v1, v2, v3, section, mat)$var
    
    ses <- c(sqrt(var_alpha), sqrt(diag(V)))
    theta <- as.vector(theta)
    sigma <- sigma
    
    list(theta = theta,
         ses = ses,
         varb = V,
         sigma = sigma,ehat = ehat, 
         effect_sch = effect_sch,
         effect_sec = effect_sec,
         effect_err = effect_err)
}



#' Linear Model Estimation using the nlme package.
#' 
#' This gets the REML or ML estimates and predictions of random effects from
#' the nlme package. 
#' function does.
#' 
#' 
#' @param x Design matrix, (p+1)xn, with intercept.  
#' @param y Response vector of nx1.
#' @param dat Data frame
#' @param method Character string indicating method to use, either "ML" or
#' "REML" (defaults to REML).
#' @return
#'   \item{theta}{ Fixed effects esimates.}
#'   \item{ses}{ Standard error for fixed effects.}
#'   \item{varb}{ Variances.}
#'   \item{sigma}{ Error. }
#'   \item{ehat}{ Raw residuals }
#'   \item{standr.lme}{ Standardized residual }
#'   \item{effect_sch}{ Cluster random error. }
#'   \item{effect_sec}{ Subcluster random error. }
#'   \item{effect_err}{ Epsilon error. }
#'   
#' @author Yusuf Bilgic
#' 
#' @seealso \code{\link{rlme}}
#' 
#' @references J. Pinheiro, D. Bates, S. DebRoy, D. Sarkar and R Development
#' Core Team. nlme linear and non- linear mixed effects models. The R Journal,
#' 2011. URL http://CRAN.R-project.org/package=nlme. R package version 3.1-98.
#' 
#' @importFrom nlme random.effects
#' @importFrom mgcv extract.lme.cov2
#' @importFrom nlme lme
#' 
#' @keywords models
#' @export
LM_est <- function(x, y, dat, method = "REML") {
    model = as.formula(paste("y ~ 1 + ", paste(colnames(x), collapse = " + ")))
    fit.lme = lme(model, data = dat, random = ~1 | school/section, method = method)
    theta0 <- extract.lme.cov2(fit.lme, dat, start.level = 3)$V[1]
    theta2 <- extract.lme.cov2(fit.lme, dat, start.level = 2)$V[[1]][1, 1] - theta0
    theta1 <- extract.lme.cov2(fit.lme, dat, start.level = 1)$V[[1]][1, 1] - (theta0 + theta2)
    sigma.l <- c(theta1, theta2, theta0)
    theta <- as.vector(summary(fit.lme)$tTable[, 1])
    intra_err.lm <- (theta0)/(sum(sigma.l))
    intra_sch.lm <- (theta1)/(sum(sigma.l))
    intra_sect.lm <- (theta1 + theta2)/(sum(sigma.l))
    ses <- as.vector(summary(fit.lme)$tTable[, 2])
    ehat <- as.vector(fit.lme$residuals[, 1])
    varb = fit.lme$varFix
    effect_sch = random.effects(fit.lme, level = 1)[, 1]
    effect_sec = random.effects(fit.lme, level = 2)[, 1]
    effect_err = as.vector(fit.lme$residuals[, 3])
    standr.lme = residuals(fit.lme, type = "pearson")
    standr.lme = as.vector(residuals(fit.lme, type = "pearson"))
    list(theta = theta,
         ses = ses,
         varb = varb,
         sigma = sigma.l, 
         ehat = ehat,
         effect_sch = effect_sch,
         effect_sec = effect_sec, 
         effect_err = effect_err,
         standr.lme = standr.lme)
}
