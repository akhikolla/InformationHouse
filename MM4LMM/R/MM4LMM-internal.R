.MM_ML <-
  function(Y , Cofactor , X , formula , Factors , VarList , Init=NULL , MaxIter = 10 , CritVar = 1e-3 , CritLogLik = 1e-3 , NbCores=1){
    
    if (is.null(Init)) Init <- rep(1,length(VarList))
    if (is.null(names(VarList))){
      NamesSigma <- paste0("Var",1:length(VarList))
    }else{
      NamesSigma <- names(VarList)
    }
    
    if (is.null(X)){
      
      
      Mat <- Cofactor
      Tmp <- model.matrix(terms(formula,keep.order=TRUE),data=as.data.frame(Mat))
      QR <- qr(Tmp)
      Fixed <- cbind(Tmp[,QR$pivot[1:QR$rank]])
      NamesFixed <- colnames(Fixed)
      res <- .MLMM(Y, Fixed , VarList , Init , MaxIter, CritVar, CritLogLik)
      
      names(res$Beta) <- NamesFixed
      rownames(res$VarBeta) <- NamesFixed
      names(res$Sigma2) <- NamesSigma
      res$attr <- attr(Tmp,"assign")[QR$pivot[1:QR$rank]]
      res$Factors <- Factors
      Res <- list(res)
      names(Res) <- "NullModel"
    }else{
      if (is.list(X)){
        Res <-  mclapply(X , function(x) {
          if (is.null(colnames(x))){
            colnames(x) <- paste0("X",1:ncol(X))
          }
          
          Mat <- cbind(Cofactor,x)
          Tmp <- model.matrix(terms(formula,keep.order=TRUE),data=as.data.frame(Mat))
          QR <- qr(Tmp)
          Fixed <- cbind(Tmp[,QR$pivot[1:QR$rank]])
          NamesFixed <- colnames(Fixed)
          
          res <- .MLMM(Y, Fixed , VarList , Init , MaxIter, CritVar, CritLogLik)
          
          names(res$Beta) <- NamesFixed
          rownames(res$VarBeta) <- NamesFixed
          names(res$Sigma2) <- NamesSigma
          res$attr <- attr(Tmp,"assign")[QR$pivot[1:QR$rank]]
          res$Factors <- Factors
          return(res)
        } ,mc.cores=NbCores)
        names(Res) <- names(X)
      }else{
        Res <-  mclapply(1:ncol(X) , function(i) {
          Mat <- Cofactor
          Mat$Xeffect <- X[,i]
          Tmp <- model.matrix(terms(formula,keep.order=TRUE),data=as.data.frame(Mat))
          #colnames(Tmp) <- gsub("Xeffect",i,colnames(Tmp))
          QR <- qr(Tmp)
          Fixed <- cbind(Tmp[,QR$pivot[1:QR$rank]])
          NamesFixed <- colnames(Fixed)
          
          res <- .MLMM(Y, Fixed , VarList , Init , MaxIter, CritVar, CritLogLik)
          
          names(res$Beta) <- NamesFixed
          rownames(res$VarBeta) <- NamesFixed
          names(res$Sigma2) <- NamesSigma
          res$attr <- attr(Tmp,"assign")[QR$pivot[1:QR$rank]]
          res$Factors <- Factors
          return(res)
        } ,mc.cores=NbCores)
        names(Res) <- colnames(X)
      }
    }
    return(Res)
  }

.MM_ML2Mat <-
  function(Y , Cofactor , X , formula , Factors , VarList , Init=NULL , MaxIter = 100 , CritVar = 1e-3 , CritLogLik = 1e-3 , NbCores=1){
    
    if (is.null(Init)) Init <- rep(1,2)
    
    K1 <- VarList[[1]]
    K2 <- VarList[[2]]
    List <- .PrepMat(Y , K1 , K2)
    Ytilde <- List$Ytilde
    Passage <- t(List$U)
    
    if (is.null(names(VarList))){
      NamesSigma <- paste0("Var",1:length(VarList))
    }else{
      NamesSigma <- names(VarList)
    }
    
    if (is.null(X)){
      Mat <- Cofactor
      Tmp <- model.matrix(terms(formula,keep.order=TRUE),data=as.data.frame(Mat))
      QR <- qr(Tmp)
      Fixed <- cbind(Tmp[,QR$pivot[1:QR$rank]])
      NamesFixed <- colnames(Fixed)
      
      res <- .MM_ML2MatRcpp(Ytilde , Fixed , Passage , List$Diag , Init ,MaxIter, CritVar , CritLogLik)
      
      names(res$Beta) <- NamesFixed
      rownames(res$VarBeta) <- NamesFixed
      names(res$Sigma2) <- NamesSigma
      res$attr <- attr(Tmp,"assign")[QR$pivot[1:QR$rank]]
      res$Factors <- Factors
      Res <- list(res)
      names(Res) <- "NullModel"
    }else{
      if (is.list(X)){
        Res <-  mclapply(X , function(x) {
          if (is.null(colnames(x))){
            colnames(x) <- paste0("X",1:ncol(X))
          }
          
          Mat <- cbind(Cofactor,x)
          Tmp <- model.matrix(terms(formula,keep.order=TRUE),data=as.data.frame(Mat))
          QR <- qr(Tmp)
          Fixed <- cbind(Tmp[,QR$pivot[1:QR$rank]])
          NamesFixed <- colnames(Fixed)
          
          res <- .MM_ML2MatRcpp(Ytilde , Fixed , Passage , List$Diag , Init ,MaxIter, CritVar , CritLogLik)
          
          names(res$Beta) <- NamesFixed
          rownames(res$VarBeta) <- NamesFixed
          names(res$Sigma2) <- NamesSigma
          res$attr <- attr(Tmp,"assign")[QR$pivot[1:QR$rank]]
          res$Factors <- Factors
          return(res)
        } ,mc.cores=NbCores)
        names(Res) <- names(X)
      }else{
        Res <-  mclapply(1:ncol(X) , function(i) {
          Mat <- Cofactor
          Mat$Xeffect <- X[,i]
          Tmp <- model.matrix(terms(formula,keep.order=TRUE),data=as.data.frame(Mat))
          #colnames(Tmp) <- gsub("Xeffect",i,colnames(Tmp))
          QR <- qr(Tmp)
          Fixed <- cbind(Tmp[,QR$pivot[1:QR$rank]])
          NamesFixed <- colnames(Fixed)
          
          res <- .MM_ML2MatRcpp(Ytilde , Fixed , Passage , List$Diag , Init ,MaxIter, CritVar , CritLogLik)
          
          names(res$Beta) <- NamesFixed
          rownames(res$VarBeta) <- NamesFixed
          names(res$Sigma2) <- NamesSigma
          res$attr <- attr(Tmp,"assign")[QR$pivot[1:QR$rank]]
          res$Factors <- Factors
          return(res)
        } ,mc.cores=NbCores)
        names(Res) <- colnames(X)
      }
    }
    return(Res)
  }

.MM_Reml <-
  function(Y , Cofactor , X , formula , Factors , VarList , Init=NULL , MaxIter = 100 , CritVar = 1e-3 , CritLogLik = 1e-3 , NbCores=1){
    
    if (is.null(Init)) Init <- rep(1,length(VarList))
    if (is.null(names(VarList))){
      NamesSigma <- paste0("Var",1:length(VarList))
    }else{
      NamesSigma <- names(VarList)
    }
    if (is.null(X)){
      Mat <- Cofactor
      Tmp <- model.matrix(terms(formula,keep.order=TRUE),data=as.data.frame(Mat))
      QR <- qr(Tmp)
      Fixed <- cbind(Tmp[,QR$pivot[1:QR$rank]])
      NamesFixed <- colnames(Fixed)
      
      res <- .RemlMM(Y, Fixed , VarList , Init , MaxIter, CritVar , CritLogLik)
      
      names(res$Beta) <- NamesFixed
      rownames(res$VarBeta) <- NamesFixed
      names(res$Sigma2) <- NamesSigma
      res$attr <- attr(Tmp,"assign")[QR$pivot[1:QR$rank]]
      res$Factors <- Factors
      Res <- list(res)
      names(Res) <- "NullModel"
    }else{
      if (is.list(X)){
        Res <-  mclapply(X , function(x) {
          if (is.null(colnames(x))){
            colnames(x) <- paste0("X",1:ncol(X))
          }
          
          Mat <- cbind(Cofactor,x)
          Tmp <- model.matrix(terms(formula,keep.order=TRUE),data=as.data.frame(Mat))
          QR <- qr(Tmp)
          Fixed <- cbind(Tmp[,QR$pivot[1:QR$rank]])
          NamesFixed <- colnames(Fixed)
          
          res <- .RemlMM(Y, Fixed , VarList , Init , MaxIter, CritVar , CritLogLik)
          
          names(res$Beta) <- NamesFixed
          rownames(res$VarBeta) <- NamesFixed
          names(res$Sigma2) <- NamesSigma
          res$attr <- attr(Tmp,"assign")[QR$pivot[1:QR$rank]]
          res$Factors <- Factors
          return(res)
        } ,mc.cores=NbCores)
        names(Res) <- names(X)
      }else{
        Res <-  mclapply(1:ncol(X) , function(i) {
          Mat <- Cofactor
          Mat$Xeffect <- X[,i]
          Tmp <- model.matrix(terms(formula,keep.order=TRUE),data=as.data.frame(Mat))
          #colnames(Tmp) <- gsub("Xeffect",i,colnames(Tmp))

          QR <- qr(Tmp)
          Fixed <- cbind(Tmp[,QR$pivot[1:QR$rank]])
          NamesFixed <- colnames(Fixed)
          
          res <- .RemlMM(Y, Fixed , VarList , Init , MaxIter, CritVar , CritLogLik)
          
          names(res$Beta) <- NamesFixed
          rownames(res$VarBeta) <- NamesFixed
          names(res$Sigma2) <- NamesSigma
          res$attr <- attr(Tmp,"assign")[QR$pivot[1:QR$rank]]
          res$Factors <- Factors
          return(res)
        } ,mc.cores=NbCores)
        names(Res) <- colnames(X)
      }
    }
    return(Res)
  }

.MM_Reml1Mat <-
  function(Y , Cofactor , X , formula , Factors , VarInv , logdetVar , NamesSigma , NbCores=1){
    
    if (is.null(NamesSigma)){
      NamesSigma <- "Var"
    }
    if (is.null(X)){
      Mat <- Cofactor
      Tmp <- model.matrix(terms(formula,keep.order=TRUE),data=as.data.frame(Mat))
      QR <- qr(Tmp)
      Fixed <- cbind(Tmp[,QR$pivot[1:QR$rank]])
      NamesFixed <- colnames(Fixed)
      
      res <- .RemlMM1Mat(Y, Fixed , VarInv , logdetVar)
      
      names(res$Beta) <- NamesFixed
      rownames(res$VarBeta) <- NamesFixed
      names(res$Sigma2) <- NamesSigma
      res$attr <- attr(Tmp,"assign")[QR$pivot[1:QR$rank]]
      res$Factors <- Factors
      Res <- list(res)
      names(Res) <- "NullModel"
    }else{
      if (is.list(X)){
        Res <-  mclapply(X , function(x) {
          if (is.null(colnames(x))){
            colnames(x) <- paste0("X",1:ncol(X))
          }
          
          Mat <- cbind(Cofactor,x)
          Tmp <- model.matrix(terms(formula,keep.order=TRUE),data=as.data.frame(Mat))

          QR <- qr(Tmp)
          Fixed <- cbind(Tmp[,QR$pivot[1:QR$rank]])
          NamesFixed <- colnames(Fixed)
          
          res <- .RemlMM1Mat(Y, Fixed , VarInv , logdetVar)
          
          names(res$Beta) <- NamesFixed
          rownames(res$VarBeta) <- NamesFixed
          names(res$Sigma2) <- NamesSigma
          res$attr <- attr(Tmp,"assign")[QR$pivot[1:QR$rank]]
          res$Factors <- Factors
          return(res)
        } ,mc.cores=NbCores)
        names(Res) <- names(X)
      }else{
        Res <-  mclapply(1:ncol(X) , function(i) {
          Mat <- Cofactor
          Mat$Xeffect <- X[,i]
          Tmp <- model.matrix(terms(formula,keep.order=TRUE),data=as.data.frame(Mat))
          #colnames(Tmp) <- gsub("Xeffect",i,colnames(Tmp))

          QR <- qr(Tmp)
          Fixed <- cbind(Tmp[,QR$pivot[1:QR$rank]])
          NamesFixed <- colnames(Fixed)
          
          res <- .RemlMM1Mat(Y, Fixed , VarInv , logdetVar)
          
          names(res$Beta) <- NamesFixed
          rownames(res$VarBeta) <- NamesFixed
          names(res$Sigma2) <- NamesSigma
          res$attr <- attr(Tmp,"assign")[QR$pivot[1:QR$rank]]
          res$Factors <- Factors
          return(res)
        } ,mc.cores=NbCores)
        names(Res) <- colnames(X)
      }
    }
    return(Res)
  }


.MM_Reml2Mat <-
  function(Y , Cofactor , X , formula , Factors , VarList , Init=NULL , MaxIter = 100 , CritVar = 1e-3 , CritLogLik = 1e-3 , NbCores=1){
    
    if (is.null(Init)) Init <- rep(1,1)
    
    K1 <- VarList[[1]]
    K2 <- VarList[[2]]
    List <- .PrepMat(Y , K1 , K2)
    Ytilde <- List$Ytilde
    Passage <- t(List$U)
    
    
    
    
    if (is.null(names(VarList))){
      NamesSigma <- paste0("Var",1:length(VarList))
    }else{
      NamesSigma <- names(VarList)
    }
    if (is.null(X)){
      Mat <- Cofactor
      Tmp <- model.matrix(terms(formula,keep.order=TRUE),data=as.data.frame(Mat))
      QR <- qr(Tmp)
      Fixed <- cbind(Tmp[,QR$pivot[1:QR$rank]])
      NamesFixed <- colnames(Fixed)
      
      res <- .MM_Reml2MatRcpp(Ytilde , Fixed , Passage , c(List$Diag) , Init ,MaxIter, CritVar , CritLogLik)
      
      names(res$Beta) <- NamesFixed
      rownames(res$VarBeta) <- NamesFixed
      names(res$Sigma2) <- NamesSigma
      res$attr <- attr(Tmp,"assign")[QR$pivot[1:QR$rank]]
      res$Factors <- Factors
      Res <- list(res)
      names(Res) <- "NullModel"
    }else{
      if (is.list(X)){
        Res <-  mclapply(X , function(x) {
          if (is.null(colnames(x))){
            colnames(x) <- paste0("X",1:ncol(X))
          }
          
          Mat <- cbind(Cofactor,x)
          Tmp <- model.matrix(terms(formula,keep.order=TRUE),data=as.data.frame(Mat))

          QR <- qr(Tmp)
          Fixed <- cbind(Tmp[,QR$pivot[1:QR$rank]])
          NamesFixed <- colnames(Fixed)
          
          res <- .MM_Reml2MatRcpp(Ytilde , Fixed , Passage , c(List$Diag) , Init ,MaxIter, CritVar , CritLogLik)
          
          names(res$Beta) <- NamesFixed
          rownames(res$VarBeta) <- NamesFixed
          names(res$Sigma2) <- NamesSigma
          res$attr <- attr(Tmp,"assign")[QR$pivot[1:QR$rank]]          
          res$Factors <- Factors
          return(res)
        } ,mc.cores=NbCores)
        names(Res) <- names(X)
      }else{
        Res <-  mclapply(1:ncol(X) , function(i) {
          Mat <- Cofactor
          Mat$Xeffect <- X[,i]
          Tmp <- model.matrix(terms(formula,keep.order=TRUE),data=as.data.frame(Mat))
          #colnames(Tmp) <- gsub("Xeffect",i,colnames(Tmp))

          QR <- qr(Tmp)
          Fixed <- cbind(Tmp[,QR$pivot[1:QR$rank]])
          NamesFixed <- colnames(Fixed)
          
          res <- .MM_Reml2MatRcpp(Ytilde , Fixed , Passage , c(List$Diag) , Init ,MaxIter, CritVar , CritLogLik)
          
          names(res$Beta) <- NamesFixed
          rownames(res$VarBeta) <- NamesFixed
          names(res$Sigma2) <- NamesSigma
          res$attr <- attr(Tmp,"assign")[QR$pivot[1:QR$rank]]      
          res$Factors <- Factors
          return(res)
        } ,mc.cores=NbCores)
        names(Res) <- colnames(X)
      }
    }
    return(Res)
  }

.MM_RemlHen <-
  function(Y , Cofactor , X , formula , Factors , Z , GList , GinvList , Rinv , logdetV , NameVar=NULL , Init=NULL , MaxIter = 100 , CritVar = 1e-3 , CritLogLik = 1e-3 , NbCores=1){
    
    if (is.null(NameVar)){
      NamesSigma <- paste0("Var",1:length(NameVar))
    }else{
      NamesSigma <- NameVar
    }
    if (is.null(X)){
      Mat <- Cofactor
      Tmp <- model.matrix(terms(formula,keep.order=TRUE),data=as.data.frame(Mat))
      QR <- qr(Tmp)
      Fixed <- cbind(Tmp[,QR$pivot[1:QR$rank]])
      NamesFixed <- colnames(Fixed)
      
      res <- .RemlMMHen(Y, Fixed , Z , GList , GinvList , Rinv , logdetV , Init , MaxIter, CritVar , CritLogLik)
      
      names(res$Beta) <- NamesFixed
      rownames(res$VarBeta) <- NamesFixed
      names(res$Sigma2) <- NamesSigma
      res$attr <- attr(Tmp,"assign")[QR$pivot[1:QR$rank]]      
      res$Factors <- Factors
      Res <- list(res)
      names(Res) <- "NullModel"
    }else{
      if (is.list(X)){
        Res <-  mclapply(X , function(x) {
          if (is.null(colnames(x))){
            colnames(x) <- paste0("X",1:ncol(X))
          }
          
          Mat <- cbind(Cofactor,x)
          Tmp <- model.matrix(terms(formula,keep.order=TRUE),data=as.data.frame(Mat))

          QR <- qr(Tmp)
          Fixed <- cbind(Tmp[,QR$pivot[1:QR$rank]])
          NamesFixed <- colnames(Fixed)
          
          res <- .RemlMMHen(Y, Fixed , Z , GList , GinvList , Rinv , logdetV , Init , MaxIter, CritVar , CritLogLik)
          
          names(res$Beta) <- NamesFixed
          rownames(res$VarBeta) <- NamesFixed
          names(res$Sigma2) <- NamesSigma
          res$attr <- attr(Tmp,"assign")[QR$pivot[1:QR$rank]]    
          res$Factors <- Factors
          return(res)
        } ,mc.cores=NbCores)
        names(Res) <- names(X)
      }else{
        Res <-  mclapply(1:ncol(X) , function(i) {
          Mat <- Cofactor
          Mat$Xeffect <- X[,i]
          Tmp <- model.matrix(terms(formula,keep.order=TRUE),data=as.data.frame(Mat))
          #colnames(Tmp) <- gsub("Xeffect",i,colnames(Tmp))

          QR <- qr(Tmp)
          Fixed <- cbind(Tmp[,QR$pivot[1:QR$rank]])
          NamesFixed <- colnames(Fixed)
          
          res <- .RemlMMHen(Y, Fixed , Z , GList , GinvList , Rinv , logdetV , Init , MaxIter, CritVar , CritLogLik)
          
          names(res$Beta) <- NamesFixed
          rownames(res$VarBeta) <- NamesFixed
          names(res$Sigma2) <- NamesSigma
          res$attr <- attr(Tmp,"assign")[QR$pivot[1:QR$rank]]    
          res$Factors <- Factors
          return(res)
        } ,mc.cores=NbCores)
        names(Res) <- colnames(X)
      }
    }
    return(Res)
  }

.MM_RemlHenDiag <-
  function(Y , Cofactor , X , formula , Factors , Z , GList , GinvList , Rinv , logdetV , NameVar=NULL , Init=NULL , MaxIter = 10 , CritVar = 1e-3 , CritLogLik = 1e-3 , NbCores=1){
    
    if (is.null(NameVar)){
      NamesSigma <- paste0("Var",1:length(NameVar))
    }else{
      NamesSigma <- NameVar
    }
    if (is.null(X)){
      Mat <- Cofactor
      Tmp <- model.matrix(terms(formula,keep.order=TRUE),data=as.data.frame(Mat))
      QR <- qr(Tmp)
      Fixed <- cbind(Tmp[,QR$pivot[1:QR$rank]])
      NamesFixed <- colnames(Fixed)
      
      res <- .RemlMMHenDiag(Y, Fixed , Z , GList , GinvList , Rinv , logdetV , Init , MaxIter, CritVar , CritLogLik)
      
      names(res$Beta) <- NamesFixed
      rownames(res$VarBeta) <- NamesFixed
      names(res$Sigma2) <- NamesSigma
      res$attr <- attr(Tmp,"assign")[QR$pivot[1:QR$rank]]      
      res$Factors <- Factors
      Res <- list(res)
      names(Res) <- "NullModel"
    }else{
      if (is.list(X)){
        Res <-  mclapply(X , function(x) {
          if (is.null(colnames(x))){
            colnames(x) <- paste0("X",1:ncol(X))
          }
          
          Mat <- cbind(Cofactor,x)
          Tmp <- model.matrix(terms(formula,keep.order=TRUE),data=as.data.frame(Mat))

          QR <- qr(Tmp)
          Fixed <- cbind(Tmp[,QR$pivot[1:QR$rank]])
          NamesFixed <- colnames(Fixed)
          
          res <- .RemlMMHenDiag(Y, Fixed , Z , GList , GinvList , Rinv , logdetV , Init , MaxIter, CritVar , CritLogLik)
          
          names(res$Beta) <- NamesFixed
          rownames(res$VarBeta) <- NamesFixed
          names(res$Sigma2) <- NamesSigma
          res$attr <- attr(Tmp,"assign")[QR$pivot[1:QR$rank]]      
          res$Factors <- Factors
          return(res)
        } ,mc.cores=NbCores)
        names(Res) <- names(X)
      }else{
        Res <-  mclapply(1:ncol(X) , function(i) {
          Mat <- Cofactor
          Mat$Xeffect <- X[,i]
          Tmp <- model.matrix(terms(formula,keep.order=TRUE),data=as.data.frame(Mat))
          #colnames(Tmp) <- gsub("Xeffect",i,colnames(Tmp))

          QR <- qr(Tmp)
          Fixed <- cbind(Tmp[,QR$pivot[1:QR$rank]])
          NamesFixed <- colnames(Fixed)
          
          res <- .RemlMMHenDiag(Y, Fixed , Z , GList , GinvList , Rinv , logdetV , Init , MaxIter, CritVar , CritLogLik)
          
          names(res$Beta) <- NamesFixed
          rownames(res$VarBeta) <- NamesFixed
          names(res$Sigma2) <- NamesSigma
          res$attr <- attr(Tmp,"assign")[QR$pivot[1:QR$rank]]      
          res$Factors <- Factors
          return(res)
        } ,mc.cores=NbCores)
        names(Res) <- colnames(X)
      }
    }
    return(Res)
  }

.AnovaTest_KR <- function(ResMMEst , TestedCombination=NULL , Cofactor = NULL , X = NULL , formula = NULL , VarList , NbCores=1){
  
  if (is.null(VarList)) stop('VarList should be provided when using the KR option')

  if (is.null(Cofactor)){
    Cofactor <- matrix(rep(1,nrow(VarList[[1]])),ncol=1)
  }  
  Cofactor <- as.data.frame(Cofactor)
  
  if (!is.data.frame(Cofactor)){
    if (!is.matrix(Cofactor)){
      stop('Cofactor should be either a matrix or a data.frame')        
    } else {
      Cofactor <- as.data.frame(Cofactor)
      names(Cofactor) <- colnames(Cofactor)
    }
    Names <- names(Cofactor)
    if (is.null(Names)){
      names(Cofactor) <- paste0("Cof",1:ncol(Cofactor))
    } else {
      names(Cofactor) <- gsub('([[:punct:]])|\\s+','',Names)
    }
  }
  CofName <- names(Cofactor)
  
  
  if (is.matrix(X)){
    X <- lapply(1:ncol(X), function(i) {
      res <- matrix(X[,i],nrow=nrow(X))
      return(res)
    })
  }
  
  if (!is.null(X)){
    if (is.matrix(X)){
      if (is.null(colnames(X))) colnames(X) <- paste0("X",1:ncol(X))
      Xname <- "Xeffect"
    }else{
      if (is.list(X)){
        if (is.null(colnames(X[[1]]))){
          Xname <- paste0("X",1:ncol(X[[1]]))
        }else{
          Xname <- colnames(X[[1]])
        }
      }else{
        stop("X should be either a matrix or a list.")
      }
    }
  }else{
    Xname <- NULL
  }
  
  if (is.null(formula)){
    formulaCof <- as.formula(paste0("~",paste0(CofName,collapse="+")))
    formulaComp <- as.formula(paste0("~",paste0(c(CofName,Xname),collapse="+")))
    Factors <- c(CofName,Xname)
  }else{
    SplitForm <- strsplit(as.character(formula)," \\+ ")
    Factors <- SplitForm[[length( SplitForm)]]
    WhatInX <- unique(unlist(sapply(Xname,function(n) grep(n,Factors))))
    if (length(WhatInX)!=0){
      formulaCof <- as.formula(paste0("~",Factors[-WhatInX],collapse="+"))
    }else{
      formulaCof <- as.formula(paste0("~",Factors,collapse="+"))
    }
    formulaComp <- formula
  }
  
  
  
  
  if (!is.null(TestedCombination)){
    #message("Wald tests about the combination are computed")
    if (is.matrix(TestedCombination)) TestedCombination <- list(TestedCombination)
    if (!is.list(TestedCombination)) stop("TestedCombination should be either a list of matrices or a matrix")
    
    TestedCombination <- lapply(1:length(TestedCombination) , function(ind) {
      c <- TestedCombination[[ind]]
      QR <- qr(c)
      if (QR$rank != nrow(c)){
        message(paste0("Contrast matrix number " ,ind," was not full rank and has been reduced."))
        return(matrix(c[QR$pivot[1:QR$rank],],ncol=ncol(c)))
      }else{
        return(c)
      }
    })
    
    
    NbVar <- length(VarList)
    
    
    Res <- mclapply(1:length(ResMMEst) , function(s) {
      x <- X[[s]]
      if (is.null(colnames(x))){
        colnames(x) <- paste0("X",1:ncol(x))
      }
      Mat <- cbind(Cofactor,x)
      Tmp <- model.matrix(terms(formulaComp,keep.order=TRUE),data=as.data.frame(Mat))
      QR <- qr(Tmp)
      Fixed <- cbind(Tmp[,QR$pivot[1:QR$rank]])
      NamesFixed <- colnames(Fixed)
      
      ResMM_x <- ResMMEst[[s]]
      
      SigmaEst <- ResMM_x$Sigma2
      BETA <- ResMM_x$Beta
      Phi <- ResMM_x$VarBeta
      Sigma <- Reduce('+',lapply(1:NbVar , function(k) VarList[[k]]*SigmaEst[k]))
      Sigma_Inv <- solve(Sigma)
      VarInf_Cof <- crossprod(Sigma_Inv,Fixed)
      List_P <- lapply(1:NbVar , function(k) -crossprod(VarInf_Cof , crossprod(VarList[[k]] , VarInf_Cof)))
      List_Q <- lapply(combn(NbVar,2,simplify=F) , function(x) crossprod(VarInf_Cof , crossprod(VarList[[x[1]]],crossprod(Sigma_Inv,crossprod(VarList[[x[2]]],VarInf_Cof)))))
      names(List_Q) <- sapply(combn(NbVar,2,simplify=F) , function(x) paste0(x[1],x[2]))
      for (i in 1:NbVar){
        List_Q[[paste0(i,i)]] <- crossprod(VarInf_Cof , crossprod(VarList[[i]],crossprod(Sigma_Inv,crossprod(VarList[[i]],VarInf_Cof))))
      }
      
      SigVar <- lapply(1:NbVar , function(k) crossprod(Sigma_Inv,VarList[[k]]))
      Phi_P <- lapply(1:NbVar , function(k) crossprod(Phi,List_P[[k]]))
      W_inv <- matrix(0,NbVar,NbVar)
      for (i in 1:NbVar){
        for (j in i:NbVar){
          W_inv[i,j] <- sum(diag(crossprod(t(SigVar[[i]]),SigVar[[j]])))/2 - sum(diag(2*Phi %*% List_Q[[paste0(i,j)]] - crossprod(t(Phi_P[[i]]),Phi_P[[j]])))/2
          
        }
        for (j in 1:i){
          W_inv[i,j] <- sum(diag(crossprod(t(SigVar[[i]]),SigVar[[j]])))/2 - sum(diag(2*Phi %*% t(List_Q[[paste0(j,i)]]) - crossprod(t(Phi_P[[j]]),Phi_P[[i]])))/2
        }
      }
      W <- ginv(W_inv)
     
      Test <- sapply(TestedCombination , function(L) {
        THETA <- crossprod(L ,  solve(tcrossprod(tcrossprod(L,Phi),L) , L))
        l <- qr(L)$rank
        p <- ncol(L)
        
        Delta <- matrix(0,NbVar,NbVar)
        A1 <- 0
        A2 <- 0
        for (i in 1:NbVar){
          for (j in 1:NbVar){
            if (i < j) {
              Delta <- Delta + W[i,j] * (List_Q[[paste0(i,j)]] - List_P[[i]] %*% Phi_P[[j]])
            }else{
              Delta <- Delta + W[i,j] * (t(List_Q[[paste0(j,i)]]) - List_P[[i]] %*% Phi_P[[j]])
            }
            A1 <- A1 + W[i,j] * sum(diag(THETA %*% Phi_P[[i]] %*% Phi)) * sum(diag(THETA %*% Phi_P[[j]] %*% Phi))
            A2 <- A2 + W[i,j] * sum(diag(THETA %*% Phi_P[[i]] %*% Phi %*% THETA %*% Phi_P[[j]] %*% Phi))
          }
        }
        DELTA <- Phi %*% Delta %*% Phi
        Num <- Phi + 2 * DELTA
        g <- ((l+1)*A1 - (l+4)*A2)/((l+2)*A2)
        c1 <- g/(3*l+2*(1-g))
        c2 <- (l-g)/(3*l+2*(1-g))
        c3 <- (l+2-g)/(3*l+2*(1-g))
        
        B <- 1/(2*l) * (A1 + 6*A2)
        E <- 1/(1-A2/l)
        V <- 2/l * ((1+c1*B)/((1-c2*B)**2*(1-c3*B)))
        
        rho <- V/(2*E**2)
        m <- 4 + (l+2)/(l*rho-1)
        lambda <- m/(E*(m-2))
        
        Stat <- 1/l * t(L %*% BETA) %*% solve(L%*% Num %*% t(L) , L%*%BETA)
        Stat_Modif <- lambda * Stat
        
        pval <- 1-pf(Stat_Modif,l,m)
        
        res <- c(Stat_Modif,pval,m)
        names(res) <- c("Stat","pval","ddf")
        return(res)
      }) 
      return(Test)  
    },mc.cores=NbCores)
    names(Res) <- names(ResMMEst)
  }else{
    Res <- mclapply(1:length(ResMMEst) , function(s) {
      
      x <- X[[s]]
      if (is.null(colnames(x))){
        colnames(x) <- paste0("X",1:ncol(x))
      }
      Mat <- cbind(Cofactor,x)
      Tmp <- model.matrix(terms(formulaComp,keep.order=TRUE),data=as.data.frame(Mat))
      QR <- qr(Tmp)
      Fixed <- cbind(Tmp[,QR$pivot[1:QR$rank]])
      NamesFixed <- colnames(Fixed)
      
      ResMM_x <- ResMMEst[[s]]
      BETA <- ResMM_x$Beta
      Attr <- ResMM_x$attr
      CommonAttr <- lapply(unique(Attr) , function(y) which(Attr==y))
      names(CommonAttr) <- c("(Intercept)",Factors[unique(Attr)[-1]])
      MatTIII <- lapply(names(CommonAttr) , function(y){
        EffectSelected <- CommonAttr[[y]]
        mat <- t(sapply(EffectSelected , function(z) as.numeric(1:length(BETA)==z)))
        return(mat)				
      })
      
      SigmaEst <- ResMM_x$Sigma2
      Phi <- ResMM_x$VarBeta
      Sigma <- Reduce('+',lapply(1:NbVar , function(k) VarList[[k]]*SigmaEst[k]))
      Sigma_Inv <- solve(Sigma)
      VarInf_Cof <- crossprod(Sigma_Inv,Fixed)
      List_P <- lapply(1:NbVar , function(k) -crossprod(VarInf_Cof , crossprod(VarList[[k]] , VarInf_Cof)))
      List_Q <- lapply(combn(NbVar,2,simplify=F) , function(x) crossprod(VarInf_Cof , crossprod(VarList[[x[1]]],crossprod(Sigma_Inv,crossprod(VarList[[x[2]]],VarInf_Cof)))))
      names(List_Q) <- sapply(combn(NbVar,2,simplify=F) , function(x) paste0(x[1],x[2]))
      for (i in 1:NbVar){
        List_Q[[paste0(i,i)]] <- crossprod(VarInf_Cof , crossprod(VarList[[i]],crossprod(Sigma_Inv,crossprod(VarList[[i]],VarInf_Cof))))
      }
      SigVar <- lapply(1:NbVar , function(k) crossprod(Sigma_Inv,VarList[[k]]))
      Phi_P <- lapply(1:NbVar , function(k) crossprod(Phi,List_P[[k]]))
      W_inv <- matrix(0,NbVar,NbVar)
      for (i in 1:NbVar){
        for (j in i:NbVar){
          W_inv[i,j] <- sum(diag(crossprod(t(SigVar[[i]]),SigVar[[j]])))/2 - sum(diag(2*Phi %*% List_Q[[paste0(i,j)]] - crossprod(t(Phi_P[[i]]),Phi_P[[j]])))/2
          
        }
        for (j in 1:i){
          W_inv[i,j] <- sum(diag(crossprod(t(SigVar[[i]]),SigVar[[j]])))/2 - sum(diag(2*Phi %*% t(List_Q[[paste0(j,i)]]) - crossprod(t(Phi_P[[j]]),Phi_P[[i]])))/2
        }
      }
      W <- ginv(W_inv)
      
      Test <- sapply(MatTIII , function(L) {
        THETA <- crossprod(L ,  solve(tcrossprod(tcrossprod(L,Phi),L) , L))
        l <- qr(L)$rank
        p <- ncol(L)
        
        Theta_Phi_P_Phi <- lapply(1:NbVar , function(k) crossprod(THETA , tcrossprod(Phi_P[[k]] , Phi)))
        
        Delta <- matrix(0,NbVar,NbVar)
        A1 <- 0
        A2 <- 0
        for (i in 1:NbVar){
          for (j in 1:NbVar){
            if (i < j) {
              Delta <- Delta + W[i,j] * (List_Q[[paste0(i,j)]] - List_P[[i]] %*% Phi_P[[j]])
            }else{
              Delta <- Delta + W[i,j] * (t(List_Q[[paste0(j,i)]]) - List_P[[i]] %*% Phi_P[[j]])
            }
            A1 <- A1 + W[i,j] * sum(diag(Theta_Phi_P_Phi[[i]])) * sum(diag(Theta_Phi_P_Phi[[j]]))
            A2 <- A2 + W[i,j] * sum(diag(Theta_Phi_P_Phi[[i]] %*% Theta_Phi_P_Phi[[j]]))
          }
        }
        DELTA <- Phi %*% Delta %*% Phi
        Num <- Phi + 2 * DELTA
        g <- ((l+1)*A1 - (l+4)*A2)/((l+2)*A2)
        c1 <- g/(3*l+2*(1-g))
        c2 <- (l-g)/(3*l+2*(1-g))
        c3 <- (l+2-g)/(3*l+2*(1-g))
        
        B <- 1/(2*l) * (A1 + 6*A2)
        E <- 1/(1-A2/l)
        V <- 2/l * ((1+c1*B)/((1-c2*B)**2*(1-c3*B)))
        
        rho <- V/(2*E**2)
        m <- 4 + (l+2)/(l*rho-1)
        lambda <- m/(E*(m-2))
        
        Stat <- 1/l * t(L %*% BETA) %*% solve(L%*% Num %*% t(L)) %*% (L%*%BETA)
        Stat_Modif <- lambda * Stat
        
        pval <- 1-pf(Stat_Modif,l,m)
        
        res <- c(Stat_Modif,m,pval)
        names(res) <- c("Stat","pval","ddf")
        return(res)
      })
      return(Test)  
    },mc.cores=NbCores)
    names(Res) <- names(ResMMEst)
  }
  
  return(Res)
}

