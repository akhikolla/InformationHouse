UmxBestMatchesFromIsland <- function(Umatrix, BestMatches, Imx, Cls= NULL, Toroid=T , RemoveOcean = T){
# V <- UmxBestMatchesFromIsland(Umx, BestMatches, Imx, Cls)
# Cuts out the Imx (Island) from the Umatrix, BestMatches and Cls
  # INPUT
  # Umatrix
  # BestMatches
  # Imx
  # Optional
  # Cls
  # Toroid         should the umatrix be Toroid? (was it calculated on a toroid)
  # RemoveOcean   should the umatrix shortened to the actual island?
  # OUTPUT
  # list of
  # Umatrix
  # BestMatches
  # Cls

  # auf vierfache groesse kacheln
  if(Toroid){
    tU <- ToroidUmatrix(Umatrix, BestMatches, Cls)
    Umatrix <- tU$Umatrix
    BestMatches <- tU$BestMatches
    Cls <- tU$Cls
  }

  BestMatches = CheckBestMatches(BestMatches, Cls)
  
  # configure filter, so that every bestmatch stays in
  if(!is.null(BestMatches)){
    BestMatchesFilter = rep(T,nrow(BestMatches)) # every Bestmatch stays
  }

  # put Imx on Umatrix and bestmatches if given
  if(!is.null(Imx)){
    for(i in 1:nrow(Imx)){
      for(j in 1:ncol(Imx)){
        if(Imx[i,j] == 1){
          Umatrix[i,j] = NA
          if(!is.null(BestMatches))
            BestMatchesFilter[(BestMatches[,2] == i) & (BestMatches[,3] == j)] = F
        }
      }
    }

    if(!is.null(BestMatches)) BestMatches = BestMatches[BestMatchesFilter,]
    if((!is.null(Cls)) & (!is.null(BestMatches))) Cls = Cls[BestMatchesFilter]
  }

  #### remove ocean around Umatrix
  if(RemoveOcean){
    oceanLine = !apply(Umatrix, 1, function(x) any(x != -1))
    startLine = min(which(!oceanLine),na.rm=T)
    endLine = length(oceanLine) - min(which(rev(!oceanLine)),na.rm=T) + 1

    oceanCol = !apply(Umatrix, 2, function(x) any(x != -1))
    startCol = min(which(!oceanCol),na.rm=T)
    endCol = length(oceanCol) - min(which(rev(!oceanCol)),na.rm=T) + 1

    if(!is.null(BestMatches)){
      BestMatches <- BestMatches - cbind(rep(0,nrow(BestMatches)),startLine-1,startCol-1)
    }
    Umatrix <- Umatrix[startLine:endLine,startCol:endCol]
  }

  Umatrix[which(is.na(Umatrix))] = 0

  return(list(Umatrix = Umatrix, BestMatches = BestMatches, Cls=Cls))
}
