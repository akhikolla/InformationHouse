#' Compare Fits
#' 
#' Compares two model fits. It returns tdbeta value and cfits values of two
#' fits. The function uses the fixed effects estimates from fit 1 and fit 2
#' along with the covariance of the rank-based fit.
#' 
#' 
#' @param x Matrix of covariates
#' @param fit1 A class of type rlme.
#' @param fit2 A class of type rlme.
#' @return Returns tdbeta and cfits values.
#' @seealso \code{\link{fitdvcov}}
#' @keywords models
#' @examples
#' 
#' 
#' data(schools)
#' model = y ~ 1 + sex + age + (1 | region) + (1 | region:school)
#' 
#' # Extract covariants into matrix
#' cov = as.matrix(data.frame(schools[,"sex"], schools[,"age"]))
#' 
#' # Fit the models using each method
#' reml.fit = rlme(model, schools, method="reml")
#' gr.fit = rlme(model, schools, method="gr")
#' 
#' compare.fits(cov, reml.fit, gr.fit)
#' 
#' @export
#' 
compare.fits <- function(x, fit1, fit2) {
    cov = as.matrix(data.frame(y = rep(1, length(x[, 1])), x))
    fit1.beta = fit1$fixed.effects$Estimate
    fit2.beta = fit2$fixed.effects$Estimate
    var.b = fit1$var.b
    tdbeta = fitdvcov(cov, fit1.beta, fit2.beta, var.b)$tdbeta
    cfits = fitdvcov(cov, fit1.beta, fit2.beta, var.b)$cfits
    return(tdbeta)
    return(cfits)
}
