createMaskFromPolygon <- function(Height, Width, CutoutPol){
  # Sets every value of the Mask to 0, that is not within CutoutPol or else to 1
  # INPUT
  # Height                            Height of umatrix
  # Width                             Width of umatrix
  # CutoutPol                         dataframe with x and y
  # OUTPUT
  # UmatrixMask[1:lines, 1:cols]      Umatrix Mask
  # author: FL

  #library(fields)
  UmatrixMask = matrix(1,nrow=Height, ncol=Width)

  # cut out polygon out of umatrix
  for(i in 1:ncol(UmatrixMask)){
    umatrixRow = cbind(rep(i,nrow(UmatrixMask)), 1:nrow(UmatrixMask))
    umatrixRowFilter = in.poly(umatrixRow, CutoutPol)
    UmatrixMask[(umatrixRowFilter), i] = 0  # set everything to 1 that is within the polygon
  }

  UmatrixMask
}
