checkUmatrixIntegrity = function(Lrn = NULL, Umx = NULL, BestMatches = NULL, Cls = NULL, Key = NULL, Imx = NULL, Wts = NULL,
                                 UStarMatrix = NULL, PMatrix = NULL){
# Checks if all given parameters are compatible to each other (Lrn and BestMatches have the same amount of rows...)
# BestMatches are expected to be given in full form. (with keys, 3 columns)

  error = c()

  tryCatch({
  ### Check DataTypes:
  print("Check DataTypes")
  if(!is.null(Lrn)) if(!is.matrix(Lrn)) return("Lrn is not a matrix")
  if(!is.null(Umx)) if(!is.matrix(Umx)) return("Umx is not a matrix")
  if(!is.null(BestMatches)) if(!is.matrix(BestMatches)) return("BestMatches is not a vector")
  if(!is.null(BestMatches)) if(ncol(BestMatches)!=3) return(paste("BestMatches have",ncol(BestMatches),"instead of 3 columns"))
  if(!is.null(Cls)) if(!is.vector(Cls)) return("Cls is not a vector")
  if(!is.null(Key)) if(!is.vector(Key)) return("Key is not a vector")
  if(!is.null(Imx)) if(!is.matrix(Imx)) return("Imx is not a matrix")
  if(!is.null(Wts)) if(!is.matrix(Wts)) return("Wts is not a matrix")
  if(!is.null(UStarMatrix)) if(!is.matrix(UStarMatrix)) return("UStarMatrix is not a matrix")
  if(!is.null(PMatrix)) if(!is.matrix(PMatrix)) return("PMatrix is not a matrix")

  ### Check LRN
  print("Check LRN")
  if(!is.null(Lrn)){
    if(!is.null(BestMatches)){
      if(nrow(Lrn) != nrow(BestMatches))
        error = c(error, "Rows of Lrn do not match the rows of BestMatches")
    }
    if(!is.null(Cls)){
      if(nrow(Lrn) != length(Cls))
        error = c(error, "Rows of Lrn do not match the rows of Cls")
    }
    if(!is.null(Key)){
      if(nrow(Lrn) != length(Key))
        error = c(error, "Rows of Lrn do not match the rows of Key")
    }
  }

  ### Check UMX
  print("Check UMX")
  if(!is.null(Umx)){
    if(!is.null(UStarMatrix)){
      if(any(dim(UStarMatrix)!=dim(Umx)))
        error = c(error, "UStarMatrix does not fit onto the Umatrix")
    }
    if(!is.null(PMatrix)){
      if(any(dim(Umx)!=dim(PMatrix)))
        error = c(error, "Umatrix does not fit onto the PMatrix")
    }
    if(!is.null(BestMatches)){
      if((nrow(Umx) < max(BestMatches[,2])) ||
         (ncol(Umx) < max(BestMatches[,3]))){
        error = c(error, "BestMatches are out of the Umatrix")
      }
    }
    if(!is.null(Imx)){
      if(any((dim(Umx)*2)!=dim(Imx)))
        error = c(error, "Imx does not fit onto the Umatrix")
    }
    if(!is.null(Wts)){
      if((nrow(Umx) * ncol(Umx)) != nrow(Wts))
        error = c(error, "Amount of neurons in wts do not match up the the size of the Umatrix")
    }
  }

  ### Check UStarMatrix
  print("Check UStarMatrix")
  if(!is.null(UStarMatrix)){
    if(!is.null(BestMatches)){
      if((nrow(UStarMatrix) < max(BestMatches[,2])) ||
         (ncol(UStarMatrix) < max(BestMatches[,3]))){
        error = c(error, "BestMatches are out of the UStarMatrix")
      }
    }
    if(!is.null(PMatrix)){
      if(any(dim(UStarMatrix)!=dim(PMatrix)))
        error = c(error, "UStarMatrix does not fit onto the PMatrix")
    }
    if(!is.null(Imx)){
      if(any((dim(UStarMatrix)*2)!=dim(Imx)))
        error = c(error, "Imx does not fit onto the UStarMatrix")
    }
    if(!is.null(Wts)){
      if((nrow(UStarMatrix) * ncol(UStarMatrix)) != nrow(Wts))
        error = c(error, "Amount of neurons in wts do not match up the the size of the UStarmatrix")
    }
  }

  ### Check PMatrix
  print("Check PMatrix")
  if(!is.null(PMatrix)){
    if(!is.null(BestMatches)){
      if((nrow(PMatrix) < max(BestMatches[,2])) ||
         (ncol(PMatrix) < max(BestMatches[,3]))){
        error = c(error, "BestMatches are out of the PMatrix")
      }
    }
    if(!is.null(Imx)){
      if(any((dim(PMatrix)*2)!=dim(Imx)))
        error = c(error, "Imx does not fit onto the PMatrix")
    }
    if(!is.null(Wts)){
      if((nrow(PMatrix) * ncol(PMatrix)) != nrow(Wts))
        error = c(error, "Amount of neurons in wts do not match up the the size of the UStarmatrix")
    }
  }

  ### Check BestMatches
  print("Check BestMatches")
  if(!is.null(BestMatches)){
    if(!is.null(Cls)){
      if(nrow(BestMatches) != length(Cls))
        error = c(error, "Number of BestMatches and of entries in Cls do not match up.")
    }
    if(!is.null(Key)){
      if(nrow(BestMatches) != length(Key))
        error = c(error, "Rows of Key do not match the rows of BestMatches")
    }
  }

  ### Check Cls
  print("Check Cls")
  if(!is.null(Cls)){
    if(!is.null(Key)){
      if(length(Cls) != length(Key))
      error = c(error, "Length of Cls do not match the rows of Key")
    }
  }

  ### Check Imx
  print("Check Imx")
  if(!is.null(Imx)){
    if(!is.null(Wts)){
      if((nrow(Imx)*ncol(Imx)/2) != nrow(Wts))
      error = c(error, "Size of Imx and amount of neurons in Wts do not match up")
    }

  }
  }, error=function(x) error = c(error, x))

  return(error)
}
