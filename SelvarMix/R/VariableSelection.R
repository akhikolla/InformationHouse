VariableSelection<-
  function(data,
           nbcluster,
           models,
           criterion,
           OrderVariable,
           hsize,
           supervised,
           z,
           nbcores)
  {
    data <- as.matrix(data)
    nbcluster <- as.integer(nbcluster)
    criterion <- as.character(criterion)
    
    
    listModels.size <- length(models["listModels"])
    nbcluster.size <- length(nbcluster)
    OutputVector.size <- listModels.size*nbcluster.size
    #print(c("listModels.size ...", listModels.size)) 
    wrapper.selectVar <- function(arg)
    {
      myModel <-  mixmodGaussianModel(listModels=models["listModels"][arg[2]])
      orderVar <- rep(NA, ncol(data))
      
      if(length(nbcluster)==1)
        orderVar <- OrderVariable 
      else
        orderVar <- OrderVariable[arg[1],]
      
      my.list <- try(rcppSelectS(data, 
                                 orderVar, 
                                 nbcluster[arg[1]], 
                                 myModel, 
                                 hsize, 
                                 criterion,
                                 z,
                                 supervised))
      if(class(my.list)=="try-error")
        return(my.list)
      else
      {
        ResSelectVar <- list()
        ResSelectVar$S <- my.list$S
        ResSelectVar$nbcluster <- my.list$nbcluster
        ResSelectVar$criterion <- my.list$criterion
        ResSelectVar$model <- my.list$model
        ResSelectVar$criterionValue  <- my.list$criterionValue
        ResSelectVar$parameters  <- my.list$parameters
        ResSelectVar$partition <- my.list$partition
        ResSelectVar$proba <- my.list$proba
        OrderAux <- setdiff(OrderVariable, ResSelectVar$S)
        ResSelectVar$W <- rcppSelectW(data, OrderAux, ResSelectVar$S, hsize)
        return(ResSelectVar)
      }
    }
    
    if(OutputVector.size < nbcores)
      nbcores <- OutputVector.size
    
    arg.grid <- matrix(0, OutputVector.size, 2)  
    arg.grid <- as.matrix(expand.grid(1:nbcluster.size, 1:listModels.size))
    ## si on est sous windows
    # junk <- list()
    # for(r in 1:nrow(arg.grid))
    # { 
    #   junk[[r]] <- try(wrapper.selectVar(arg.grid[r,]))
    #   print(paste0("r = ", r))
    #   print(junk[[r]])
    #   if(class(junk[[r]])=="try-error")
    #     print(junk[[r]])
    # }  
    if(Sys.info()["sysname"] == "Windows")
    {
      cl <- makeCluster(nbcores)
      common.objects <- c("data", 
                          "OrderVariable", 
                          "nbcluster",
                          "models",
                          "hsize", 
                          "criterion",
                          "supervised",
                          "z")
      clusterExport(cl=cl, varlist = common.objects, envir = environment())
      clusterEvalQ(cl, require(Rmixmod))
      junk <- parApply(cl = cl,  
                       X = arg.grid,
                       MARGIN = 1,
                       FUN = wrapper.selectVar)
      stopCluster(cl)
      
    }
    else
      junk <- mclapply(X = as.list(data.frame(t(arg.grid))),
                       FUN = wrapper.selectVar,
                       mc.cores = nbcores,
                       mc.silent = FALSE,
                       mc.preschedule = TRUE,
                       mc.cleanup = TRUE)
    nb.fails <- 0
    for(idx in 1:OutputVector.size)
      if(class(junk[[idx]]) == "try-error")
        nb.fails <- nb.fails + 1
    
    VariableSelectRes <-  vector(length = (OutputVector.size - nb.fails), mode ="list")
    
    idx <- 1
    for(ll in 1:OutputVector.size)
      if(class(junk[[ll]])!="try-error")
      {
        VariableSelectRes[[idx]]$S <- junk[[ll]][["S"]]
        VariableSelectRes[[idx]]$W <- junk[[ll]][["W"]]
        VariableSelectRes[[idx]]$U <- setdiff(1:dim(data)[2], union(junk[[ll]][["S"]], junk[[ll]][["W"]]))
        VariableSelectRes[[idx]]$criterionValue <- junk[[ll]][["criterionValue"]]
        VariableSelectRes[[idx]]$criterion <- junk[[ll]][["criterion"]]
        VariableSelectRes[[idx]]$model <- junk[[ll]][["model"]]
        VariableSelectRes[[idx]]$nbcluster <- junk[[ll]][["nbcluster"]]
        VariableSelectRes[[idx]]$parameters <- junk[[ll]][["parameters"]]
        VariableSelectRes[[idx]]$partition <- junk[[ll]][["partition"]]
        VariableSelectRes[[idx]]$proba <- junk[[ll]][["proba"]]
        idx <- idx + 1
      }
    
    return(VariableSelectRes)
    
  }


