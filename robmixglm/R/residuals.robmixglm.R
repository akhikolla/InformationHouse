residuals.robmixglm <-
  function(object,
           type = c("deviance", "pearson"),
           ...)
  {
    deviance.residuals <- function(){
      calc.binomial <- function(){
        n <- object$Y[,1]+object$Y[,2]
        sign(object$Y[,1]-n*fitted(object))*sqrt(
          2*(object$Y[,1]*log(object$Y[,1]/(n*fitted(object)))+
             object$Y[,2]*log(object$Y[,2]/(n*(1-fitted(object))))))
      }
      res <- switch(object$family,
            gaussian=(object$Y-fitted(object))/sqrt(coef(object$fit)[length(coef(object$fit))]),
            binomial=calc.binomial(),
            poisson=sign(object$Y-fitted(object))*sqrt(2*(object$Y*log(object$Y/fitted(object))-(object$Y-fitted(object)))),
            gamma=sign(object$Y-fitted(object))*sqrt(2*(-log(object$Y/fitted(object))+(object$Y-fitted(object))/fitted(object))),
            truncpoisson=sign(object$Y-fitted(object))*sqrt(2*(object$Y*log(object$Y/fitted(object))-(object$Y-fitted(object))))
      )
      res
    }
    pearson.residuals <- function(){
      res <- switch(object$family,
           gaussian=(object$Y-fitted(object))/sqrt(coef(object$fit)[length(coef(object$fit))]),
           binomial=(object$Y[,1]-(object$Y[,1]+object$Y[,2])*fitted(object))/
             sqrt((object$Y[,1]+object$Y[,2])*fitted(object)*(1-fitted(object))),
           poisson=(object$Y-fitted(object))/sqrt(fitted(object)),
           gamma=(object$Y-fitted(object))/fitted(object),
           truncpoisson=sign(object$Y-fitted(object))/(sqrt(fitted(object)*exp(fitted(object))*(-1+exp(fitted(object))-fitted(object)))/(exp(fitted(object))-1))
      )
      res
    }
    type <- match.arg(type)
     res <- switch(type,
           deviance=deviance.residuals(),
           pearson=pearson.residuals()
     )
    res
  }