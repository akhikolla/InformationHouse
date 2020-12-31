CheckInputsL <- function(x, z, lambda, rho, type, hsize, models, rmodel, imodel, xtest, ztest, nbcores)
{
  if(missing(x)){
    stop("x is missing !")
  } 
  if(is.matrix(x) == FALSE & is.data.frame(x) == FALSE) 
    stop(paste(sQuote("x"), "must be a matrix"))

  if(missing(lambda)){
    stop("lambda is missing!")
  }
  if(is.vector(lambda) == FALSE | length(lambda) <= 1){ 
    stop(paste(sQuote("lambda"), "must be a vector with length >= 2"))
  }
  if (sum(lambda<=0)){
    stop("lambda must be greater than 0!")
  }
  
  if(missing(rho)){
    stop("rho is missing!")
  }
  if(is.vector(rho) == FALSE){ 
    stop(paste(sQuote("rho"), "must be a vector"))
  }
  if(sum(rho<=0)){
    stop("rho must be greater than 0!")
  }
  
  if(missing(hsize)){
    hsize <- 3
  }
  if(!is.wholenumber(hsize) | sum(hsize < 1) | hsize > ncol(x)) 
    stop(paste(sQuote("hsize"), "must be a positive integer <= ncol(x)"))
  
  if(missing(models)){
    models <- mixmodGaussianModel(listModels = c("Gaussian_pk_L_C", 
                                                 "Gaussian_pk_Lk_C", 
                                                 "Gaussian_pk_L_Ck", 
                                                 "Gaussian_pk_Lk_Ck"))
  }
  if(!(isS4(models) && is(models, "GaussianModel")))
  {
    stop("models must be an GaussianModel S4 object! (see Rmixmod package)")
  }
  
  if(missing(rmodel)){
    rmodel <- c("LI", "LB", "LC")
  }
  if( sum(rmodel %in% c("LI","LB","LC")) != length(rmodel) ){
    stop(cat(rmodel[which(!(rmodel %in% c("LI","LB","LC")))], "is not a valid rmodel name !\n"))
  }
  
  if(missing(imodel)){
    imodel <- c("LI", "LB")
  }
  if ( sum(imodel %in% c("LI","LB")) != length(imodel) ){
    stop(cat(imodel[which(!(imodel %in% c("LI","LB")))], "is not a valid imodel name !\n"))
  }
  
  if( missing(z)){
    stop("labels are missing!")
  }
  
  if(is.factor(z))
  {
    z <- factor(z, labels = seq(1,length(levels(z))))
    z <- as.integer(z)
  }
  
  if (min(z) <= 0 | length(z) != nrow(x)){
    stop("Each observation in z must have a valid cluster affectation !")
  }

  if(missing(xtest) && !missing(ztest))
    stop("xtest is missing!")
  
  if(!missing(xtest) && missing(ztest))
    stop("ztest are missing!")
  
  if(!missing(xtest))
    if((is.matrix(xtest) == FALSE & is.data.frame(xtest) == FALSE)) 
      stop(paste(sQuote("xtest"), "must be a matrix"))
  
 
  if(missing(xtest) && missing(ztest))
    testing <- FALSE
  
  if((is.wholenumber(nbcores) == FALSE) || (nbcores < 1)) 
    stop(paste(sQuote("nbcores"), "must be an integer > 0"))
}  