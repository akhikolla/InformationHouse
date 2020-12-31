UmatrixProportions <- function(Umatrix, Imx, Toroid){
  if(is.null(Umatrix)) return(list(Width = 80, Height = 50))

  if(!Toroid){
    Width = ncol(Umatrix)
    Height = nrow(Umatrix)
  }
  else if(!is.null(Imx)){
    lineHeights = rowSums(Imx,na.rm=T)
    islandTop = min(which(lineHeights<ncol(Imx)),na.rm=T)
    islandBottom = length(lineHeights) - min(which(rev(lineHeights)<ncol(Imx)),na.rm=T) + 1

    colHeights = colSums(Imx,na.rm=T)
    islandLeft = min(which(colHeights<nrow(Imx)),na.rm=T)
    islandRight = length(colHeights) - min(which(rev(colHeights)<nrow(Imx)),na.rm=T) + 1

    Width = islandRight - islandLeft + 1
    Height = islandBottom - islandTop + 1
  }
  else{
    Width = ncol(Umatrix)*2
    Height = nrow(Umatrix)*2
  }

  return(list(Width = Width, Height = Height))
}
