SortvarLearn <- function(x,
                         z,
                         type="lasso",
                         lambda=seq(20, 100, by = 10),
                         rho=seq(1, 2, length=2),
                         nbcores=min(2,  detectCores(all.tests = FALSE, logical = FALSE)))
{
  # check x parameter
  if(missing(x)){
    stop("x is missing !")
  } 
  if(is.matrix(x) == FALSE && is.data.frame(x) == FALSE) 
    stop(paste(sQuote("x"), "must be a matrix"))  
  
  # check lambda parameter
  #if(missing(lambda)){
  #  stop("lambda is missing!")
  #}
  if(is.vector(lambda) == FALSE || length(lambda) <= 1){ 
    stop(paste(sQuote("lambda"), "must be a vector with length >= 2"))
  }
  if (sum(lambda<=0)){
    stop("lambda must be greater than 0!")
  }
  
  
  # check rho parameter
  #if(missing(rho)){
  #  stop("rho is missing!")
  #}
  if(is.vector(rho) == FALSE){ 
    stop(paste(sQuote("rho"), "must be a vector"))
  }
  if(sum(rho<=0)){
    stop("rho must be greater than 0!")
  }
  
  
  # check whether the z is missing
  if ( missing(z)){
    stop("labels are missing!")
  }
  
  
  if(missing(z) || length(z)==0){
    warning("z are missing, intrumental initialization without output effect")
    z <- rep(1, nrow(x)) 
    
  }
  
  
  if(min(z) <= 0 || length(z) != nrow(x)){
    stop("Each observation in z must have a valid cluster affectation !")
  }
  
  # check nbcores 
  if((is.wholenumber(nbcores) == FALSE) || (nbcores < 1)) 
    stop(paste(sQuote("nbcores"), "must be an integer > 0"))
  
  x <- as.matrix(scale(x, TRUE, TRUE))
  n <- as.integer(nrow(x))
  p <- as.integer(ncol(x))
  nbcluster <- as.integer(max(z))
  if(type=="lasso")
  {
  VarRole <- matrix(NA,(length(lambda)*length(rho)), p) 
  VarRole <- DiscriminantAnalysisGlasso(x, nbcluster, lambda, rho, knownlabels = z, nbcores)
  var.role.sum <- colSums(VarRole) 
  OrderVariable <- sort.int(var.role.sum,decreasing=TRUE,index.return=TRUE)$ix
  }
  if(type=="likelihood")
      OrderVariable <-orderlikL(x, z, nbcores)
  
  return(OrderVariable)    
}
