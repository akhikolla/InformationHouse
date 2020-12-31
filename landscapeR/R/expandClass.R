#' Expand an existing class of patches.
#
#' @inheritParams makePatch
#' @inheritParams makeClass
#' @param class The raster value of class (or patch) to expand.
#' @param size integer. Size of expansion, as number of raster cells.
#' @param bgr integer. The background available where expansion is allowed (i.e. shrinking classes).
#' @return A RasterLayer object. If \code{rast=FALSE} returns a list of vectors, each containing the \code{context} cells assigned to each patch.
#' @examples
#' library(raster)
#'
#' m = matrix(0, 33, 33)
#' r = raster(m, xmn=0, xmx=10, ymn=0, ymx=10)
#' r = makeClass(r, 5, 10)
#' plot(r)
#'
#' rr = expandClass(r, 1, 100)
#' plot(rr)
#'
#' ## This function can be used to mimic shapes, by providing a skeleton:
#' m[,17] = 1
#' r = raster(m, xmn=0, xmx=10, ymn=0, ymx=10)
#' plot(r)
#'
#' rr = expandClass(r, 1, 100)
#' plot(rr)
#' @export
expandClass <- function(context, class, size, bgr=0, pts = NULL) {
  if(length(class) > 1){ stop('A single value only is admitted for "class" argument.') }
  if(class %in% bgr){ warning('Value to attribute to patches same as background cells value (arg. "class" equals "bgr").') }
  if(any(is.na(size) | size <=0)){ stop('Invalid "size" argument provided.') }
  bd <- raster::boundaries(context, type='outer', classes=TRUE, directions=8)
  bd <- t(raster::as.matrix(bd))
  if(!is.matrix(context)) {
    mtx <- t(raster::as.matrix(context))
  } else {
    mtx <- context
  }
  if(is.null(pts)){
    edg <- which(bd==1 & mtx==class)
  } else {
    edg <- pts
  }

  if(length(bgr) > 1){
    p_bgr <- which(is.element(mtx, bgr[-1]))
    vals <- mtx[p_bgr]
    bgr <- bgr[1]
    bgrCells <- c(which(mtx == bgr), p_bgr)
    .assignValues(bgr, p_bgr, mtx) # mtx[p_bgr] <- bgr
  } else {
    bgrCells <- which(mtx == bgr)
  }
  if(length(bgrCells) == 0){stop('No cells available, landscape full')}
  if(size > (length(bgrCells))){stop('Expansion size bigger than available landscape')}
  if(length(edg) == 1){
    pts <- edg
  } else {
    pts <- sample(edg, 1)
  }
  cg <- 1
  while(cg < size){
    ad <- .contigCells(pts, bgr, mtx)
    if(length(ad) == 0) {
      edg <- edg[edg != pts]
      if(length(edg) <= 1) {
        if(cg == 1){
          warning('Expanding classes do not touch shrinking classes. Input raster returned.')
          break
        } else {
          warning('Maximum patch size reached: ', cg)
          break
        }
      }
      pts <- sample(edg, 1)
      next
    }
    .assignValues(class, ad, mtx) # mtx[ad] <- class
    cg <- cg + length(ad)
    edg <- c(edg[edg != pts], ad)
    if(length(edg) == 1){
      pts <- edg
    } else {
      pts <- sample(edg, 1)
    }
  }
  if(exists('p_bgr')){
#    id <- mtx[p_bgr] != class
    id <- .subset(mtx, p_bgr) != class
#    mtx[p_bgr[id]] <- vals[id]
    mtx[.subset(p_bgr, id)] <- .subset(vals, id)
  }
  context[] <- t(mtx)
  return(context)
}
