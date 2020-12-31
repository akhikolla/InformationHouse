orderlikC<- function(x, g, nbcores)
{
  
  loglikwrapper <-function(j,x,g){return(mixmodCluster(x[,j], nbCluster=g, models=mixmodGaussianModel(family = "all"))@bestResult@likelihood)}
  
  return(order(simplify2array(mclapply(X=1:ncol(x), 
                                       FUN = loglikwrapper, 
                                       x=x,
                                       nbcluster=g,
                                       mc.preschedule = TRUE, 
                                       mc.silent = TRUE, 
                                       mc.cleanup = TRUE, 
                                       mc.cores = min(nbcores, ncol(x)))), decreasing = TRUE))
}