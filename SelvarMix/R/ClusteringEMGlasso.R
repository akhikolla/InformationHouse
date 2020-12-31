ClusteringEMGlasso <- function(data, 
                               nbcluster, 
                               lambda, 
                               rho,
                               nbcores)
{
  data <- as.matrix(data)
  n <- as.integer(dim(data)[1])
  p <- as.integer(dim(data)[2])
  nbcluster <- as.integer(nbcluster)
  
  
  if((length(lambda)*length(rho)) < nbcores)
    nbcores <- (length(lambda)*length(rho))
  
  if(Sys.info()["sysname"] == "Windows")
    cl <- makeCluster(nbcores)
  
  if(length(nbcluster) == 1)
    
    junk <- InitParameter(data, nbcluster, n.start = 250, small.pen = 0.5) 
  
  else{  
    wrapper.init.parameter <- function(k){return(InitParameter(data, k, n.start = 250, small.pen = 0.5))}
    if(Sys.info()["sysname"] == "Windows")
    {
      common.objects <- c("InitParameter")
      clusterEvalQ(cl, require(glasso))
      clusterExport(cl=cl, varlist = common.objects, envir = environment())
      junk <- clusterApply(cl, x = as.integer(nbcluster), fun = wrapper.init.parameter)
    }
    else
      junk <- mclapply(X = as.integer(nbcluster), 
                       FUN = wrapper.init.parameter, 
                       mc.cores = nbcores,
                       mc.preschedule = TRUE,
                       mc.cleanup = TRUE)
    
    
  }
  
  
  wrapper.clusteringEMGlasso <- function(prm)
  {
    result <- rcppClusteringEMGlasso(P, prm[1], prm[2])
    return(result)
  }
  ##...................................................................##
  ##...................................................................##
  pen.grid <- matrix(0, (length(lambda)*length(rho)), 2)  
  pen.grid <- as.matrix(expand.grid(lambda, rho))
  VarRole <- array(0,dim=c((length(lambda)*length(rho)), p, length(nbcluster)))
  parallel.varrole <- list()
  if(length(nbcluster)==1)
  {
    P <- junk
    if(Sys.info()["sysname"] == "Windows")
    {
      common.objects <- c("P") 
      clusterEvalQ(cl, require(glasso))
      clusterExport(cl=cl, varlist = common.objects, envir = environment())
      parallel.varrole[[1]] <-  parApply(cl, 
                                          X = pen.grid,
                                          MARGIN =  1,
                                          FUN = wrapper.clusteringEMGlasso)  
      
    }
    else
      parallel.varrole[[1]] <-  mclapply(X = as.list(data.frame(t(pen.grid))), 
                                         FUN = wrapper.clusteringEMGlasso,
                                         mc.cores = nbcores,
                                         mc.preschedule = TRUE,
                                         mc.cleanup = TRUE)
  }
  else
    for(k in 1:length(nbcluster))
    {
      P <- junk[[k]]
      if(Sys.info()["sysname"] == "Windows")
      {
        common.objects <- c("P")
        clusterEvalQ(cl, require(glasso))
        clusterExport(cl=cl, varlist = common.objects, envir = environment())
        parallel.varrole[[k]] <- parApply(cl, 
                                          X = pen.grid,
                                          MARGIN =  1,
                                          FUN = wrapper.clusteringEMGlasso)
      }
      else
        parallel.varrole[[k]] <- mclapply(X = as.list(data.frame(t(pen.grid))),
                                          FUN = wrapper.clusteringEMGlasso,
                                          mc.cores = nbcores,
                                          mc.preschedule = TRUE,
                                          mc.cleanup = TRUE)
    } 
  
  if(Sys.info()["sysname"] == "Windows")
    stopCluster(cl)
  for(k in 1:length(nbcluster))
  {
    var.role <- matrix(0,(length(lambda)*length(rho)), p)
    for(j in 1:nrow(var.role))
      if(class(parallel.varrole[[k]][[j]])!="try-error")
        var.role[j,] <- parallel.varrole[[k]][[j]]   
    
    VarRole[,,k] <- var.role
  }
  
  return(VarRole)
}







