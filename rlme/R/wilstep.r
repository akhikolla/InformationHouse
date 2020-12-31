#' Wilcoxon One Step Rank-based Estimate in GR Method
#' 
#' Gets weighted rank based fittings for nested designs.
#' 
#' Initial inputs are from the independent model.
#' 
#' @param I Number of clusters.
#' @param sec A vector of subcluster numbers in clusters.
#' @param mat A matrix of numbers of observations in subclusters.  Dimension is
#' Ixmax(number ofsubclusters). Each row indicates one cluster.
#' @param init boolean
#' @param y Response vector of nx1.
#' @param x Design matrix, pxn, without intercept.
#' @param sigmaa2 Initial sigma for cluster in three-level design. 
#' @param sigmaw2 Initial sigma for subcluster in three-level design.
#' @param sigmae2 Initial sigma for error in three-level design.
#' @param thetaold Initial input. 
#' @param eps Epsilon value 
#' @param iflag2 y or n 
#' @param rprpair Either 'hl-disp' or 'med-mad'
#' 
#' @author J. W. McKean and Y. K. Bilgic
#' 
#' @references Y. K. Bilgic and J. W. McKean. Iteratively reweighted
#' generalized rank-based method in mixed models. 2013. Under preperation.
#' 
#' J. T. Terpstra and J. W. McKean. Rank-based analysis of linear models using
#' R. Journal of Statistical Software, 14(7) 1 - 26, 7 2005. ISSN 1548-7660.
#' URL http://www.jstatsoft.org/v14/i07.
#' 
#' @export
wilstep <- function(I, sec, mat, init = F, y, x, sigmaa2 = 1, sigmaw2 = 1, 
    sigmae2 = 1, thetaold = c(0), eps = 1e-04, iflag2 = 0, rprpair = "hl-disp") {
    location = scale = 2
    if (rprpair == "med-mad") {
        location = scale = 1
    }
    n = length(y)
    if (init == T) {
        #fitw = wwest(x, y, print.tbl = F)
        fit = minimize_dispersion(as.matrix(x), as.matrix(y))
      
        #theta = fitw$tmp1$coef
        #ehat = fitw$tmp1$residuals
        
        theta = fit$theta
        ehat = fit$ehat
        
        fitvc = rprmeddis(I, sec, mat, ehat, location, scale, 
            rprpair = rprpair)
        sigmaa2 = fitvc$siga2
        sigmaw2 = fitvc$sigw2
        sigmae2 = fitvc$sigmae2
        rea = fitvc$frei
        rew = fitvc$frew
        ree = fitvc$free
    }
    else {
        sigmay = sigymake(I, sec, mat, sigmaa2, sigmaw2, sigmae2)
        sigma12inv = sigmay$sigy12i
        ystar = sigma12inv %*% y
        xstar = sigma12inv %*% cbind(rep(1, n), x)
        
        #fitw = wwest(xstar, ystar, print.tbl = F)
        
        fit = minimize_dispersion(as.matrix(xstar), as.matrix(ystar))
        
        # Subtract intercept from residuals
        fit.residuals = fit$ehat - fit$theta[1]
        
        #yhat = ystar - fitw$tmp1$residuals
        yhat = ystar - fit.residuals
        
        fitcsp = projcsp(xstar, yhat)
        theta = fitcsp$betahat
        ehat = y - cbind(rep(1, n), x) %*% theta
        fitvc = rprmeddis(I, sec, mat, ehat, location, scale, 
            rprpair = rprpair)
        sigmaa2 = fitvc$siga2
        sigmaw2 = fitvc$sigw2
        sigmae2 = fitvc$sigmae2
        rea = fitvc$frei
        rew = fitvc$frew
        ree = fitvc$free
        chk = sum((theta - thetaold)^2)/sum(thetaold^2)
        if (chk < eps) {
            iflag2 = 1
        }
    }
    list(theta = theta, ehat = ehat, sigmaa2 = sigmaa2, sigmaw2 = sigmaw2, 
        sigmae2 = sigmae2, rea = rea, rew = rew, ree = ree, iflag2 = iflag2)
}
