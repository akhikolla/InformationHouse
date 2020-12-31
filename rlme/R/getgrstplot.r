#' Q-Q Plot and Standardized Residual Plot for the GR fit.
#' 
#' It gets Q-Q Plot and Standardized Residual Plot of residuals.
#' 
#' The fit is obtained from rlme()
#' 
#' 
#' @param rlme.fit RLME fit object
#' @seealso rlme 
#' 
#' @export
getgrstplot <- function(rlme.fit) {
    y = rlme.fit$y
    xstar = rlme.fit$xstar
    ystar = rlme.fit$ystar
    residstar.gr = rlme.fit$ehats
    ystarhat.gr = (ystar - residstar.gr)
    yhat.gr = y - rlme.fit$effect.err
    temp <- stanresidgr(x = xstar, y = ystar, resid = residstar.gr, 
        delta = 0.8, param = 2, conf = 0.95)
    standr.gr <- temp$stanr
    trim = 2
    par(mfrow = c(1, 2), font.main = 1)
    plot(standr.gr ~ yhat.gr, pch = "o", xlab = "Fit", ylab = "Standardized Residual", 
        main = "Stand. Residuals vs. Fits in GR")
    abline(h = c(-trim, trim), col = "red")
    qqnorm(standr.gr, pch = "o", main = "Normal Q-Q Plot")
    qqline(standr.gr, col = "red")
    par(mfrow = c(1, 1))
    list(sresid=standr.gr)
}



#' Q-Q Plot and Standardized Residual Plot for the REML or ML fit.
#' 
#' It gets Q-Q Plot and Standardized Residual Plot of residuals.
#' concise (1-5 lines) description of what the function does.
#' 
#' 
#' @param rlme.fit The fit is obtained from rlme()
#' 
#' @seealso rlme 
#' 
#' @export
#' 
getlmestplot <- function(rlme.fit) {
    standr.lme = rlme.fit$standr.lme
    y = rlme.fit$y
    yhat.lme = y - rlme.fit$effect.err
    trim = 2
    par(mfrow = c(1, 2), font.main = 1)
    plot(standr.lme ~ yhat.lme, pch = "o", xlab = "Fit", ylab = "Standardized Residual", 
        main = "Stand. Residuals vs. Fits in REML")
    abline(h = c(-trim, trim), col = "red")
    qqnorm(standr.lme, pch = "o", main = "Normal Q-Q Plot")
    qqline(standr.lme, col = "red")
    par(mfrow = c(1, 1))
    list(sresid=standr.lme)
}