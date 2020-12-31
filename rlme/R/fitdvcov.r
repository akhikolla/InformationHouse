#' Fitdvcov
#' 
#' Obtains measurement for the fits based on estimates beta1, beta2 and
#' covariance matrix from a rank based methods.
#' 
#' 
#' @param x1 data
#' @param beta1 model 1 beta estimate
#' @param beta2 model 2 beta estimate
#' @param vcw variance matrix
#' @seealso \code{\link{compare.fits}}
#' @examples
#' 
#' 
#' # Compare GR and JR methods
#' 
#' data(schools)
#' 
#' model = y ~ 1 + sex + age + (1 | region) + (1 | region:school)
#' 
#' # Extract covariants into matrix
#' cov = as.matrix(data.frame(schools[,"sex"], schools[,"age"]))
#' 
#' # Fit the models using each method
#' jr.fit = rlme(model, schools, method="jr")
#' gr.fit = rlme(model, schools, method="gr")
#' 
#' # Extract beta estimates, ignoring the intercept
#' jr.beta = jr.fit$fixed.effects$Estimate[c(2, 3)]
#' gr.beta = gr.fit$fixed.effects$Estimate[c(2, 3)]
#' 
#' # Extract beta variance matrix
#' var.b = jr.fit$var.b
#' 
#' fitdvcov(cov, jr.beta, gr.beta, var.b)
#' 
#' @export
fitdvcov <- function(x1, beta1, beta2, vcw) {
    n = dim(x1)[1]
    p = dim(x1)[2]
    bd = beta1 - beta2
    tdbeta = t(bd) %*% solve(vcw) %*% bd
    bmtd = (4 * (p + 1)^2)/n
    fit1 = x1 %*% beta1
    fit2 = x1 %*% beta2
    xv = x1 %*% vcw %*% t(x1)
    cfits = rep(0, n)
    for (i in 1:n) {
        cfits[i] = (fit1[i] - fit2[i])/sqrt(xv[i, i])
    }
    bmcf = 2 * sqrt((p + 1)/n)
    list(tdbeta = c(tdbeta), bmtd = bmtd, cfits = c(cfits), bmcf = bmcf)
}
