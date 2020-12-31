#' HBR Weight
#' 
#' Calculates hbr weights for the GEER method. This turns a vector of weights
#' for a vector of errors. Used to make factor space more robust, up to 50\%
#' breakdown. See HM (2012) and Terpstra and McKean (2005) for details. The ww
#' package produces this weights as well.
#' 
#' The ww package explains how it is obtained. 
#' 
#' @param xmat Design matrix, pxn, without intercept.
#' @param y Response vector in nx1.
#' @param percent This is 0.95.
#' @param intest This is obtained from myltsreg(xmat, y)$coef 
#' 
#' @author J. W. McKean
#' 
#' @seealso GEER_est
#' 
#' @references T. P. Hettmansperger and J. W. McKean. Robust Nonparametric
#' Statistical Methods. Chapman Hall, 2012.
#' 
#' J. T. Terpstra and J. W. McKean. Rank-based analysis of linear models using
#' R. Journal of Statistical Software, 14(7):1 - 26, 7 2005. ISSN 1548-7660.
#' URL http://www.jstatsoft.org/v14/i07. 
#' 
#' @importFrom robustbase covMcd
#' @importFrom MASS ltsreg
#' 
#' @export
hbrwts_gr <- function(xmat, y, percent = 0.95, intest = ltsreg(xmat,y)$coef) {
    # modified hbrwts function in ww and rlme.v2
    # need to install two packs:
    # install.packages("robustbase")
    # install.packages("MASS")  
    # library(robustbase)
    # library(MASS)
    # xmat= x, y is dependent
    # check it with older ww function
    
    xmat = as.matrix(xmat)
    
    aa = covMcd(xmat)

    if(dim(xmat)[2]==1){  
      robdis2=mahalanobis(xmat, aa$raw.center, aa$raw.cov)
    } else {
      robdis2=aa$raw.mah
    }
    
    y = as.matrix(y)
    n = dim(xmat)[1]
    p = dim(xmat)[2]
    cut = qchisq(percent, p)
    resids = y - intest[1] - xmat %*% as.matrix(intest[2:(p + 1)])
    sigma = mad(resids)
    m = psi(cut/robdis2)
    a = resids/(sigma * m)
    c = (median(a) + 3 * mad(a))^2
    h = sqrt(c)/a
    ans = psi(abs(h))
    #    ans=psi(abs(h)^sqrt(sum(h^2)))
    return(ans)
  }
