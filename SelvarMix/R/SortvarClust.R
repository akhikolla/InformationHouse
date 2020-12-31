SortvarClust <- function(x,
                         nbcluster,
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
  
  
  # check nbCluster parameter
  if(missing(nbcluster)){
    stop("nbcluster is missing!")
  }
  if(sum(!is.wholenumber(nbcluster))){
    stop("nbcluster must contain only integer!")
  }
  if(sum(nbcluster < 1)){ 
    stop(paste(sQuote("nbcluster"), "must be an integer greater than 0!"))
  }
  
  
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
  
  if((is.wholenumber(nbcores) == FALSE) || (nbcores < 1)) 
    stop(paste(sQuote("nbcores"), "must be an integer > 0"))
  
  
  x <- scale(x)
  n <- as.integer(nrow(x))
  p <- as.integer(ncol(x))
  K <- as.integer(nbcluster)
  ## Ordre des variables
  OrderVariable <- matrix(NA, nrow=length(nbcluster),ncol=p)
  if(type=="lasso")
  {
    VarRole <- array(NA,dim=c((length(lambda)*length(rho)), p, length(nbcluster))) 
    VarRole <- ClusteringEMGlasso(x,nbcluster,lambda,rho, nbcores)
    Matrix0 <- matrix(0, nrow=length(nbcluster), ncol=p)
    for (k in 1:length(nbcluster))
      Matrix0[k,]<- colSums(VarRole[,,k])    
    for (k in 1:length(nbcluster))
      OrderVariable[k,] <- sort.int(Matrix0[k,],decreasing=TRUE,index.return=TRUE)$ix
  }
  if(type=="likelihood")
    for(k in 1:length(nbcluster))
      OrderVariable[k,] <-orderlikC(x, nbcluster[k], nbcores)
  
  return(OrderVariable)    
}
