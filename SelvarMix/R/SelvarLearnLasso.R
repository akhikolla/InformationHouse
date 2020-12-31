SelvarLearnLasso <- 
  function(x,  
           z,
           lambda = seq(20, 100, by = 10), 
           rho = seq(1, 2, length=2),
           type = "lasso",
           rank,
           hsize = 3,  
           models = mixmodGaussianModel(listModels = c("Gaussian_pk_L_C", "Gaussian_pk_Lk_C", "Gaussian_pk_L_Ck", "Gaussian_pk_Lk_Ck")),
           rmodel = c("LI", "LB", "LC"),
           imodel = c("LI", "LB"),
           xtest=x,
           ztest=z, 
           nbcores = min(2,  detectCores(all.tests = FALSE, logical = FALSE)))
  {
    options(warn=-1)
    testing <- TRUE
    CheckInputsL(x, z, lambda, rho, type, hsize, models, rmodel, imodel, xtest, ztest, nbcores)
    supervised <- TRUE
    nbcluster <- as.integer(max(z))
    criterion <- "BIC" 
    x <- as.matrix(x)
    n <- as.integer(nrow(x))
    p <- as.integer(ncol(x))
    OrderVariable <- rep(NA, p) 
    xstd <- scale(x, TRUE, TRUE)
    if(missing(rank))
    {
      cat("variable  ranking\n")
      OrderVariable <- SortvarLearn(xstd, z, type, lambda, rho, nbcores)
    }
    else
      OrderVariable <- rank
    bestModel <- list()
    cat("SRUW  selection with BIC criterion\n")
    VariableSelectRes <- VariableSelection(x,
                                           nbcluster,
                                           models,
                                           criterion,
                                           OrderVariable,
                                           hsize,
                                           supervised,
                                           z,
                                           nbcores)
    #print(VariableSelectRes[[1]])
    cat("model selection with BIC criterion\n")
    bestModel <- ModelSelectionClust(VariableSelectRes,
                                     x,
                                     rmodel,
                                     imodel,
                                     nbcores)
    
    if(testing)
    {
      xAux <- as.data.frame(x[,bestModel$S])
      xtestAux <- as.data.frame(xtest[,bestModel$S])
      model <- mixmodGaussianModel(listModels = bestModel$model)
      learn <- mixmodLearn(xAux, knownLabels = z, models = model)
      predict <- mixmodPredict(xtestAux, classificationRule = learn["bestResult"])
      bestModel$proba <- predict["proba"]
      bestModel$partition <- predict["partition"]
      bestModel$error <- 1 - mean(predict["partition"] == ztest)
    }
    else
    {
      bestModel$proba <- NULL
      bestModel$partition <- NULL
      bestModel$error <- NULL
    }
    
    
    colnames(bestModel$regparameters) = bestModel$U
    rownames(bestModel$regparameters) = c("intercept",bestModel$R)
    output <- list(S=bestModel$S,
                   R=bestModel$R,
                   U=bestModel$U,
                   W=bestModel$W,
                   criterionValue=bestModel$criterionValue,
                   criterion=bestModel$criterion,
                   model=bestModel$model,
                   rmodel=bestModel$rmodel,
                   imodel=bestModel$imodel,
                   parameters=bestModel$parameters,
                   #nbcluster=bestModel$nbcluster,
                   partition=bestModel$partition,
                   proba=bestModel$proba,
                   error=bestModel$error,
                   regparameters=bestModel$regparameters)
    class(output) <- "selvarmix"
    return(output)
  
  }