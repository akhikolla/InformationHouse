orderlikL<- function(x, z, nbcores)
{
  g <- length(unique(z))
  loglikwrapper <-function(j,x,g){return(mixmodLearn(x[,j], knownLabels=z, models=mixmodGaussianModel(family = "all"))@bestResult@likelihood)}
  return(order(simplify2array(mclapply(X=1:ncol(x), 
                                       FUN = loglikwrapper, 
                                       x=x,
                                       nbcluster=g,
                                       mc.preschedule = TRUE, 
                                       mc.silent = TRUE, 
                                       mc.cleanup = TRUE, 
                                       mc.cores = min(nbcores, ncol(x)))), decreasing = TRUE))
}