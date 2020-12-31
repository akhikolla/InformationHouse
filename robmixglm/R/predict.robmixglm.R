predict.robmixglm <-
  function(object, newdata = NULL, type = c("link", "response"), ...)
{
  if (!inherits(object, 'robmixglm'))
    stop("Primary argument much be a robmixglm object")

    call <- match.call()

    if (missing(type)) type <- "link"
    
     if(is.null(newdata)) {
       lp <- object$X %*% matrix(coef(object),ncol=1) + object$offset
       if (type=="link") return(lp)
       else return(switch(
         object$family,
         gaussian = lp,
         binomial = 1.0/(1.0+exp(-lp)),
         poisson = exp(lp),
         gamma = exp(lp),
         truncpoisson = exp(lp),
         nbinom = exp(lp)
       ))
     }
    Terms <- delete.response(terms(object))
    m <- model.frame(Terms, newdata, 
                     xlev = object$xlevels)
    X <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
    offset <- rep(0, nrow(X))
    # if (!is.null(off.num <- attr(terms(object), "offset")))
    #   for(i in off.num)
    #     offset <- offset + eval(attr(terms(object), "variables")[[i+1]], newdata)
    if (!is.null(object$call$offset))
      offset <- offset + eval(object$call$offset, newdata)
    
     lp <- as.vector(X %*% matrix(coef(object),ncol=1))+offset
     if (type=="link") return(lp)
     else return(switch(
       object$family,
       gaussian = lp,
       binomial = 1.0/(1.0+exp(-lp)),
       poisson = exp(lp),
       gamma = exp(lp),
       truncpoisson = exp(lp),
       nbinom = exp(lp)
     ))
  }
