MMEst <-
  function(Y , Cofactor=NULL , X=NULL , formula = NULL , VarList , ZList=NULL , Method="Reml" , Henderson = NULL , Init=NULL , CritVar=1e-3 , CritLogLik = 1e-3 , MaxIter=100 , NbCores=1){
    
    if (Sys.info()[['sysname']]=="Windows"){
      if (NbCores!=1){
        NbCores <- 1
        message("NbCores > 1 is not supported on Windows (mclapply), NbCores is set to 1")
      }
    }
    
    #### Check on dimension
    
    
    ## Comparison Y and Cofactor
    Nind <- length(Y)
    if (is.null(Cofactor)){
      Cofactor <- matrix(rep(1,Nind),ncol=1)
    }
    if (is.vector(Cofactor)) Cofactor <- matrix(Cofactor,ncol=1)
    if ((length(Cofactor)!=0)*(Nind != nrow(Cofactor))){
      stop("Incompatible dimension between Y and Cofactor")
    }
    
    
    ## Comparison Y and X
    if ((length(X)!=0)){
      if (is.matrix(X)){
        if (nrow(X)!=Nind){
          stop("Incompatible dimension between Y and X")
        }
      }else{
        if (class(X)=="list"){
          invisible(mclapply(X , function(x) {
            if (nrow(x)!=Nind){
              stop("Incompatible dimension between Y and a matrix in the list X")
            }
          },mc.cores=NbCores))
        }else{
          stop("X should be a matrix or a list")
        }
      }
    }

    ## Comparison Y and Var
    #First check for list length
    NbVar <- length(VarList)
    NbZ <- length(ZList)
    if ((NbVar != NbZ)*(NbZ>0)){ stop("Lists ZList and VarList should have the same number of elements")}
    
    
    #### Work on fixed effect quanti/quali:
    
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
      Factors <- SplitForm[[length(SplitForm)]]
      WhatInX <- unique(unlist(sapply(Xname,function(n) grep(n,Factors))))
      if (length(WhatInX)!=0){
        formulaCof <- as.formula(paste0("~",Factors[-WhatInX],collapse="+"))
      }else{
        formulaCof <- as.formula(paste0("~",Factors,collapse="+"))
      }
      formulaComp <- formula
    }
    
    
    ## Comparaison Init and Var
    if ((length(Init)!=0)*(length(Init)!=NbVar)){
      stop("The length of Init and the number of variance are not compatible")
    }
    
    
    ## Choose of the MM program
    
    if (NbVar==1){
	
	if (NbZ > 0){
		Rdiag <- isDiagonal(ZList[[NbZ]])&&isDiagonal(VarList[[NbZ]])
		if (Rdiag){
			Tmp <- diag(ZList[[1]])^2*diag(VarList[[1]])
        		VarInv <- 1/VarList[[NbZ]]
			logdetVar <- sum(log(Tmp))
		}else{
			Tmp <- .sym_inverseRcpp(tcrossprod( tcrossprod(ZList[[1]] , VarList[[1]]) ,ZList[[1]]))
			VarInv <- Tmp$inverse
			logdetVar <- Tmp$log_det
		}
	}else{
		Tmp <- .sym_inverseRcpp(VarList[[1]])
		VarInv <- Tmp$inverse
		logdetVar <- Tmp$log_det
	}
	rm(Tmp)
	NamesSigma <- names(VarList)
	Res <- .MM_Reml1Mat(Y=Y , Cofactor = Cofactor , X = X , formula = formulaComp , Factors=Factors , VarInv=VarInv , logdetVar=logdetVar , NamesSigma = NamesSigma , NbCores=NbCores)
    }else{

    ## Choose of algorithm
    DimVar <- sapply(1:(length(VarList)-1) , function(x) ncol(VarList[[x]]) )
    if (is.null(Henderson)){
      if (Method=="ML"){
        Henderson <- FALSE
      }else{
        Henderson <- FALSE
        if ((length(Y) > sum(DimVar) + ncol(Cofactor))) Henderson <- TRUE
      }	    	
    }
    


    if (  (NbVar==2) || (!Henderson) || (Method=="ML") ){
      
      
      ## If NbZ>0 create VarList
      if (NbZ>0){
        NamesVar <- names(VarList)
        VarList <- mclapply(1:NbZ , function(x) {
          Z <- ZList[[x]]
          Var <- VarList[[x]]
          if (ncol(Z)==nrow(Var)){
            if(nrow(Z)==Nind){
              return( tcrossprod(Z,tcrossprod(Z,Var)) )
            } else {
              stop("Incompatible dimensions between Y and ZList")
            }
          } else {
            stop("Incompatible dimensions between Zlist and VarList")
          }
        },mc.cores=NbCores)
        names(VarList) <- NamesVar
      }else{
        invisible(lapply(1:length(VarList) , function(x) {
          if (nrow(VarList[[x]])!=Nind){
            stop(paste0("Incompatible dimensions between Y and VarList"))
          }
        }))
      }
      
      
      if (Method=="Reml"){
        if (NbVar==2){
          if(Henderson) message("Henderson's trick is not implemented for 2 variance components, joint diagonalization trick is used")
          InfMethod <- .MM_Reml2Mat
        }else{
          
          InfMethod <- .MM_Reml
        }
      }else{
        if (Henderson) message("Henderson's trick is not available for ML inference, classical MM is used")
        if (NbVar==2){
          InfMethod <- .MM_ML2Mat
        }else{
          InfMethod <- .MM_ML
        }
      }
      ## Creation Init
      if (length(Init)==0){
        Init <- rep(var(Y)/NbVar,NbVar)
        if (!is.null(X)){
          Res <- InfMethod(Y=Y , Cofactor=Cofactor , X=NULL , formula = formulaCof , Factors=Factors , VarList = VarList , Init = Init , CritVar=CritVar , CritLogLik = CritLogLik , MaxIter=MaxIter , NbCores=1)
          Init <- Res[[1]]$Sigma
        }
      }
      
      
      #### Launch of the MM program
      Res <- InfMethod(Y=Y , Cofactor=Cofactor , X=X , formula = formulaComp , Factors=Factors , VarList = VarList , Init = Init , CritVar=CritVar , CritLogLik = CritLogLik , MaxIter=MaxIter , NbCores=NbCores)
      # if (!is.null(X)){
      #   names(Res) <- names(ListFixed)
      # }
    }else{
      
      invisible(sapply(1:NbZ , function(x) {
        if (ncol(ZList[[x]])!=nrow(VarList[[x]])) stop("Incompatible dimensions between Zlist and VarList")
        if (nrow(ZList[[x]])!=length(Y))  stop("Incompatible dimensions between Y and ZList")
      }))
      
      
      Rdiag <- isDiagonal(ZList[[NbZ]])&&isDiagonal(VarList[[NbZ]])
      
      if (Rdiag){
        VarList[[NbZ]] <- diag(ZList[[NbZ]])^2*diag(VarList[[NbZ]])
        Rinv <- 1/VarList[[NbZ]]
        Tmp <- mclapply(1:(NbZ-1) , function(x) .sym_inverseRcpp(VarList[[x]]),mc.cores=NbCores)
        GinvList <- lapply(Tmp,function(x) x$inverse)
        logdetV <- c(sapply(Tmp,function(x) x$log_det),sum(log(VarList[[NbZ]])))
        rm(Tmp)
        GList <- VarList
        NamesVar <- names(VarList)
        logdetV <- c(logdetV,sum(log(VarList[[NbZ]])))
        rm(VarList)
        InfMethod <- .MM_RemlHenDiag
      }else{
        VarList[[NbZ]] <- tcrossprod( tcrossprod(ZList[[NbZ]] , VarList[[NbZ]]) ,ZList[[NbZ]])
        
        Tmp <- mclapply(VarList , function(x) .sym_inverseRcpp(x) , mc.cores=NbCores)
        Varinv <- lapply(Tmp , function(x) x$inverse)
        logdetV <- sapply(Tmp , function(x) x$log_det)
        Rinv <- Varinv[[length(Varinv)]]
        GinvList <- lapply(1:(length(Varinv)-1) , function(x) Varinv[[x]])
        rm(Varinv)
        GList <- VarList
        NamesVar <- names(VarList)
        rm(VarList)
        rm(Tmp)
        InfMethod <- .MM_RemlHen
      }
      
      Zg = Reduce(cbind,lapply(1:(length(ZList)-1) , function(x) ZList[[x]]))
      
      ## Creation Init
      if (length(Init)==0){
        Init <- rep(var(Y)/NbVar,NbVar)
        if (!is.null(X)){
          Res <- InfMethod(Y=Y , Cofactor = Cofactor , X=NULL , formula = formulaCof , Factors=Factors , Z = Zg , GList=GList , GinvList = GinvList , Rinv = Rinv , logdetV = logdetV , NameVar = NamesVar , Init = Init , CritVar=CritVar , CritLogLik = CritLogLik , MaxIter=MaxIter , NbCores=1)
          Init <- Res[[1]]$Sigma2
        }
      }
      
      
      #### Launch of the MM program
      Res <- InfMethod(Y=Y , Cofactor = Cofactor , X=X , formula = formulaComp , Factors=Factors , Z = Zg , GList=GList , GinvList = GinvList , Rinv = Rinv , logdetV = logdetV , NameVar = NamesVar , Init = Init , CritVar=CritVar , CritLogLik = CritLogLik , MaxIter=MaxIter , NbCores=NbCores)
      # if (!is.null(X)){
      #   names(Res) <- names(ListFixed)
      # }
      

    }}
    
    return(Res)
  }
