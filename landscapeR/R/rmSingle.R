#' Remove single tones from patches
#'
#' @description Patch creation algorithm can occasionally leave single cells scattered within patches.
#' This function reduces the "salt-pepper" effect, identifying or correcting those cells.
#' % ADD FUNCTIONALITY NOT ONLY SINGLES BUT USER DEFINED SIZE
#' @param rst input raster landscape.
#' @param rm logical, if TRUE returns the raster without single tones cells, if FALSE a vector of numbers identifying the index of single tones cells.
#' @return A raster without single tones cells. If \code{rm=FALSE}, it returns a vector of numbers identifying the index of single tones cells.
#' The value assigned to single tone cells is picked from one of the four neighbouring cells, selected at random.
#' @examples
#' library(raster)
#' m = matrix(0, 33, 33)
#' r = raster(m, xmn=0, xmx=10, ymn=0, ymx=10)
#' patchSize = 500
#'
#' ## Make a patch and introduce a single tone cell
#' r = makePatch(r, patchSize, spt=578, rast=TRUE)
#' r[578] = 0
#' plot(r)
#'
#' ## Now remove it
#' plot( rmSingle(r) )
#'
#' ## Single tones can be identified but not removed:
#' rmSingle(r, rm = FALSE)
#' @export
rmSingle <- function(rst, rm = TRUE){
  dim1 <- dim(rst)[1]
  dim2 <- dim(rst)[2]
  singles <- vector()
  v <- vval <- raster::getValues(rst)
  v <- which(is.finite(v))
  for (pt in v){  ## Faster than sapply or vapply!
    if(pt %% dim2 == 0){
      cc <- dim2
      rr <- pt / dim2
    } else {
      rr <- trunc(pt / dim2) + 1
      cc <- pt - (rr-1) * dim2
    }
    ad <- c(rr-1,rr+1,rr,rr,cc,cc,cc-1,cc+1)
    ad[ad	 <= 0 | c(ad[1:4] > dim1, ad[5:8] > dim2)] <- NA
    ad <- ad[5:8] + (ad[1:4]-1)*dim2
    ad <- ad[is.finite(ad)]
    if(all(.subset(vval, ad) != .subset(vval, pt))){
      if(rm == TRUE){
        vval[pt] <- sample(.subset(vval, ad), 1)
      } else {
        singles <- c(singles, pt)
      }
    }
  }
  if(rm == TRUE){
    rst[] <- vval
    return(rst)
  } else {
    return(singles)
  }
}
