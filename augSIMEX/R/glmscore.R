glmscore<-function(beta,Y,data,weight,offset,family){
  if (!missing(family)){
    ### makesure family is matched
    if (is.character(family))
      family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family))
      family <- family()
    if (is.null(family$family)) {
      print(family)
      stop("'family' not recognized")
    }
  } else stop("Family is missing. ")
  
  
  ### Obtain the function in the family object
  variance <- family$variance
  linkinv <- family$linkinv
  mu.eta <- family$mu.eta
  
  return(score.glm(beta,Y=Y,DataM=data,weight,offset,linkinv=linkinv,var=variance,mueta=mu.eta))
}
