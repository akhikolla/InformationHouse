ToroidUmatrix <- function(Umatrix, BestMatches = NULL, Cls = NULL){
# Umatrix <- ToroidUmatrix(Umatrix)
# make 4 umatrices out of one. Optionally do the same for the BestMatches and their classes
# INPUT
# Umatrix(1:Lines,1:Columns)	      Umatrix to be Toroid
# OPTIONAL
# BestMatches(1:n,1:2)		      Positions of BestMatches to be plotted onto the Umatrix
# Cls(1:n)			      class identifier for the bestmatch at the given point
# OUTPUT
# list with
# Umatrix(1:(2*Lines),1:(2*Columns))
# BestMatches(1:(4*n),1:2)
# Cls(1:(4*n)))
# author: copied out of plotMatrixTopView; FL
  rows=nrow(Umatrix)
  cols=ncol(Umatrix)

  Umatrix <- Umatrix[c(1:rows,1:rows),c(1:cols,1:cols)]

  bm_keys = NULL
  if(!is.null(BestMatches)){
    # extract the keys
    if(ncol(BestMatches) == 3){
      bm_keys = BestMatches[,1]
      BestMatches = BestMatches[,c(2,3)]
    }
    else{
      bm_keys = 1:nrow(BestMatches)
    }

    bmRow <- nrow(BestMatches)
    BestMatches <- BestMatches[rep(1:bmRow,4),]
    BestMatches[(bmRow+1):(2*bmRow),1] <- BestMatches[(bmRow+1):(2*bmRow),1]+rows # unterer rechter Quadrant
    BestMatches[(2*bmRow+1):(3*bmRow),2] <- BestMatches[(2*bmRow+1):(3*bmRow),2]+cols # oberer linker Quadrant
    BestMatches[(3*bmRow+1):(4*bmRow),1] <- BestMatches[(3*bmRow+1):(4*bmRow),1]+rows # oberer rechter Quadrant
    BestMatches[(3*bmRow+1):(4*bmRow),2] <- BestMatches[(3*bmRow+1):(4*bmRow),2]+cols # oberer rechter Quadrant
  }
  # If Cls not missing, adjust
  if(!is.null(Cls)){
    Cls <- rep(Cls,4)
  }
  # reattach the keys
  if(!is.null(bm_keys))
    BestMatches = cbind(rep(bm_keys,4), BestMatches)

  list(Umatrix=Umatrix,BestMatches=BestMatches,Cls=Cls)
}
