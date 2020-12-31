SelvarClustLasso <- 
  function(x, 
           nbcluster, 
           lambda = seq(20, 100, by = 10), 
           rho=seq(1, 2, length=2),
           type="lasso",
           rank,
           hsize = 3, 
           criterion = "BIC", 
           models = mixmodGaussianModel(listModels = c("Gaussian_pk_L_C", "Gaussian_pk_Lk_C", "Gaussian_pk_L_Ck", "Gaussian_pk_Lk_Ck")),
           rmodel = c("LI", "LB", "LC"),
           imodel = c("LI", "LB"),
           nbcores = min(2,  detectCores(all.tests = FALSE, logical = FALSE)))
  {
    options(warn=-1)
    OrderVariable <- matrix(NA, length(nbcluster), ncol(x))
    CheckInputsC(x, nbcluster,lambda, rho, type, hsize, criterion, models, rmodel, imodel, nbcores)
    x <- as.matrix(x)
    n <- as.integer(nrow(x))
    p <- as.integer(ncol(x))
    nbcluster <- as.integer(nbcluster)
    OrderVariable <- matrix(NA, nrow = length(nbcluster), ncol = p) 
    xstd <- scale(x, TRUE, TRUE)
    
    supervised <- FALSE 
    knownlabels <- as.integer(1:n)
    if(missing(rank))
    {
      cat("variable  ranking\n")
      OrderVariable <- SortvarClust(xstd, nbcluster, type, lambda, rho, nbcores)
    }  
    else
      for(r in 1:nrow(OrderVariable))
        OrderVariable[r,] <- rank
    bestModel <- list()
    if(length(criterion)==1)
    {
      cat("SRUW selection with", criterion, "criterion\n")
      VariableSelectRes <- VariableSelection(x,
                                             nbcluster,
                                             models,
                                             criterion,
                                             OrderVariable,
                                             hsize,
                                             supervised,
                                             knownlabels,
                                             nbcores)
      
      if(criterion=="BIC"){
        cat("model selection  with BIC criterion\n")
        bestModel$BIC <- ModelSelectionClust(VariableSelectRes,
                                             x,
                                             rmodel,
                                             imodel,
                                             nbcores)
      }
      else
      {
        cat("model selection  with ICL criterion\n")
        bestModel$ICL <- ModelSelectionClust(VariableSelectRes,
                                             x,
                                             rmodel,
                                             imodel,
                                             nbcores)
      }
    }
    else
    {
      for(crit in criterion)
      {
        cat("SRUW selection with", crit, "criterion\n")
        VariableSelectRes <- VariableSelection(x,
                                               nbcluster,
                                               models,
                                               crit,
                                               OrderVariable,
                                               hsize,
                                               supervised,
                                               knownlabels,
                                               nbcores)
        
        cat("model selection  with", crit, " criterion\n ")
        cmd <- paste('bestModel$', crit, ' <- ModelSelectionClust(VariableSelectRes,x,rmodel,imodel,nbcores)', sep ="")
        eval(parse(text = cmd))
      }  
    }
    
    if(length(bestModel)==2)
    {
      output <- vector(mode="list",length = 2)
      for(cl in 1:2)
      {
        if(length(bestModel[[cl]]$U)!=0)
          colnames(bestModel[[cl]]$regparameters) = bestModel[[cl]]$U
        if(length(bestModel[[cl]]$R)!=0)
          rownames(bestModel[[cl]]$regparameters) = c("intercept",bestModel[[cl]]$R)
        object <- list(S=bestModel[[cl]]$S,
                       R=bestModel[[cl]]$R,
                       U=bestModel[[cl]]$U,
                       W=bestModel[[cl]]$W,
                       criterionValue=bestModel[[cl]]$criterionValue,
                       criterion=bestModel[[cl]]$criterion,
                       model=bestModel[[cl]]$model,
                       rmodel=bestModel[[cl]]$rmodel,
                       imodel=bestModel[[cl]]$imodel,
                       parameters=bestModel[[cl]]$parameters,
                       nbcluster=bestModel[[cl]]$nbcluster,
                       partition=bestModel[[cl]]$partition,
                       proba=bestModel[[cl]]$proba,
                       regparameters=bestModel[[cl]]$regparameters)
        class(object) <- "selvarmix"
        output[cl] <-object
      }
    }else
    {
      if(length(bestModel[[1]]$U)!=0)
        colnames(bestModel[[1]]$regparameters) = bestModel[[1]]$U
      if(length(bestModel[[1]]$U)!=0)
        rownames(bestModel[[1]]$regparameters) = c("intercept",bestModel[[1]]$R)
      output <- list(S=bestModel[[1]]$S,
                     R=bestModel[[1]]$R,
                     U=bestModel[[1]]$U,
                     W=bestModel[[1]]$W,
                     criterionValue=bestModel[[1]]$criterionValue,
                     criterion=bestModel[[1]]$criterion,
                     model=bestModel[[1]]$model,
                     rmodel=bestModel[[1]]$rmodel,
                     imodel=bestModel[[1]]$imodel,
                     parameters=bestModel[[1]]$parameters,
                     nbcluster=bestModel[[1]]$nbcluster,
                     partition=bestModel[[1]]$partition,
                     proba=bestModel[[1]]$proba,
                     regparameters=bestModel[[1]]$regparameters)
      class(output) <- "selvarmix"
    }
    return(output)
  }

