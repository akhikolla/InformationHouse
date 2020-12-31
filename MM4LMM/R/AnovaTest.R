AnovaTest <- 
function(ResMMEst , TestedCombination=NULL , Type = "TypeIII" , Cofactor = NULL , X = NULL , formula = NULL , VarList = NULL , NbCores=1){

	if (length(ResMMEst[[1]]) != 8) stop("ResMMEst should be the output of the function MMEst")

  if (Type == "KR"){
    return(.AnovaTest_KR(ResMMEst , TestedCombination , Cofactor , X , formula , VarList , NbCores))
  }else{
  
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
		
		Res <- mclapply(ResMMEst , function(x) {
			Beta <- x$Beta
			VarBeta <- x$VarBeta

			Test <- t(sapply(TestedCombination , function(C) {
				if (ncol(C)!=length(Beta)){
					Stat <- NA
					Pval <- NA
					df <- NA
				}else{
					rgC <- qr(C)$rank
					CBeta <- tcrossprod(C,t(Beta))
					Stat <- crossprod(CBeta, solve(tcrossprod(C,tcrossprod(C,VarBeta)), CBeta))

					Pval <- pchisq(Stat , df = rgC , lower.tail=FALSE)
					df <- rgC
				}
				return(c(Stat,Pval,df))			
			}))
			colnames(Test) <- c("Wald","pval","df")
			rownames(Test) <- names(TestedCombination)

			return(Test)
		},mc.cores=NbCores)
		names(Res) <- names(ResMMEst)	
	}else{
	  Factors <- ResMMEst[[1]]$Factors
		if (Type=="TypeI"){
			Res <- mclapply(ResMMEst , function(x) {
				Beta <- x$Beta
				Attr <- x$attr
				CommonAttr <- lapply(unique(Attr) , function(y) which(Attr==y))
				names(CommonAttr) <- c("(Intercept)",Factors[unique(Attr)[-1]])

				MatTI <- lapply(names(CommonAttr) , function(y){
					mat <- as.numeric(1:length(Beta) %in% CommonAttr[[y]])
					return(mat)				
				})

				VarBeta <- x$VarBeta
				Chol <- chol(solve(VarBeta))
				TypeIcomp <- tcrossprod(Chol,t(Beta))^2
				Test <- t(sapply(MatTI , function(y) {
					if (length(y)!=length(Beta)){
						Stat <- NA
						Pval <- NA
						df <- NA
					}else{
						Stat <- sum(y*TypeIcomp)
						Pval <- pchisq(Stat , df=sum(y) , lower.tail=FALSE)
						df = sum(y)
					}
					return(c(Stat,Pval,df))
				}))

				colnames(Test) <- c("Wald (Type I)","pval","df")
				rownames(Test) <- names(CommonAttr)
				return(Test)
			},mc.cores=NbCores)
			names(Res) <- names(ResMMEst)
		}
		if (Type=="TypeIII"){
			Res <- mclapply(ResMMEst , function(x) {
				Beta <- x$Beta
				Attr <- x$attr
				CommonAttr <- lapply(unique(Attr) , function(y) which(Attr==y))
				names(CommonAttr) <- c("(Intercept)",Factors[unique(Attr)[-1]])
				MatTIII <- lapply(names(CommonAttr) , function(y){
					EffectSelected <- CommonAttr[[y]]
					mat <- t(sapply(EffectSelected , function(z) as.numeric(1:length(Beta)==z)))
					return(mat)				
				})

				VarBeta <- x$VarBeta

				Test <- t(sapply(MatTIII , function(C) {
					if (ncol(C)!=length(Beta)){
						Stat <- NA
						Pval <- NA
						df <- NA
					}else{
						rgC <- qr(C)$rank
						CBeta <- tcrossprod(C,t(Beta))
						Stat <- crossprod(CBeta, solve(tcrossprod(C,tcrossprod(C,VarBeta)), CBeta))
						Pval <- pchisq(Stat , df = rgC , lower.tail=FALSE)
						df = rgC
					}
					return(c(Stat,Pval,df))
				}))

				colnames(Test) <- c("Wald (Type III)","pval","df")
				rownames(Test) <- names(CommonAttr)
				return(Test)
			},mc.cores=NbCores)
			names(Res) <- names(ResMMEst)
		}
		if ((Type!="TypeI")&&(Type!="TypeIII")) warning("AnovaTest computes TypeI or TypeIII test when TestedCombination is NULL")
	}
	CompNA <- sum(unlist(mclapply(Res,function(x) sum(is.na(x[,2])) , mc.cores=NbCores)))
	if (CompNA > 0) message(paste0(CompNA," tests were not performed because of incompatible dimension"))
	return(Res)
  }
}



