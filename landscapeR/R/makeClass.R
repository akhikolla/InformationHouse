#' Create a class of patches.
#'
#' @inheritParams makePatch
#' @param context Raster object, a raster of an empty landscape or a mask, indicating where the patch cannot be generated (see \code{bgr} argument).
#' @param npatch number of patches per class
#' @param size integer. The size of patches, as number of raster cells. A single integer can be provided, in which case all patches will have that size.
#' @param pts integer or matrix. The seed point location around which the patches are built (random points are given by default). It can be an integer, as indexes of the cells in the raster, or a two columns matrix indicating x and y coordinates.
#' @return A RasterLayer object, or a vector of cell numbers  if \code{rast=FALSE}.
#' @details The patches created can be contiguous, therefore resembling a single patch with size
#' equal to the sum of contiguous cells. The patches are created starting from the seed points (if provided) and iteratively sampling randomly neighbouring cells at the edge of the patch, according to von Neumann neighbourhood (four cells, aka Rook case).
#' There is a tolerance of +/- 3 cells from the patch size declared in \code{size} argument.
#' @examples
#' library(raster)
#'
#' m = matrix(0, 33, 33)
#' r = raster(m, xmn=0, xmx=10, ymn=0, ymx=10)
#' num = 5
#' size = 15
#' rr = makeClass(r, num, size)
#' plot(rr)
#'
#' ## Create a class of three patches of given size at three corners of the spatial context
#' size = c(10, 50, 200)
#' pts = c(1, 33, 1089)
#' rr = makeClass(r, 3, size, pts)
#' plot(rr)
#' @export
makeClass <- function(context, npatch, size, pts = NULL, bgr=0, edge=FALSE, rast=TRUE, val=1){
  if(length(npatch) != 1){ stop('A single integer value must be provided to argument "npatch"') }
  if(length(val) > 1){ stop('A single value only is admitted for "val" argument.') }
  if(val %in% bgr){ warning('Value to attribute to patches same as background cells value (arg. "val" equals "bgr").') }
  if(length(size) != npatch & length(size) > 1){ stop('Number of patches not matching the length of size vector') }
  if(any(is.na(size) | size <=0)){ stop('Invalid "size" argument provided.') }
  if(npatch <= 0 | npatch > length(context) | is.na(npatch)) { stop('Invalid number of patches required (e.g. more than available landscape cells). Check argument "npatch".') }
  if(rast==TRUE & edge==TRUE){
    edge=FALSE
    warning('Edge output reset to FALSE. edge=TRUE only when raster output is not required (i.e. rast=FALSE)')
  }
  mtx <- t(raster::as.matrix(context))
  if(length(size)==1){
    size <- rep(size, npatch)
  }
  bgrCells <- which(is.element(mtx, bgr))
  if(length(bgrCells) == 0){
    stop('No background cells available with value ', bgr, '. Try checking argument "bgr".')
  }
  if(length(bgr > 1)){
    bgr <- bgr[1]
    .assignValues(bgr, bgrCells, mtx) #mtx[bgrCells] <- bgr
  }
  if(is.null(pts)){
    pts <- sample(bgrCells, npatch) # bgrCells[.Internal(sample(length(bgrCells), npatch, FALSE, NULL))] #
  }
  pts <- .toCellIndex(context, pts)
  ## invalidPts <- which(pts > length(mtx) | pts < 1 | pts %% 1 != 0 | is.na(pts))
  ## if(length(invalidPts) > 0){
  if(length(pts) != npatch){ stop('Number of patches not matching number of seed points provided.') }
  invalidPts <- !is.element(pts, bgrCells)
  if(any(invalidPts)){
    if(all(invalidPts)){
      stop('All seed points invalid.')
    } else {
      warning('Invalid seed points: ', paste(pts[invalidPts], collapse='; '), '\n Invalid points were ignored.')
      pts <- pts[!invalidPts]
      size <- size[!invalidPts]
      npatch <- length(pts)
    }
  }
  lst <- list()
  for(np in 1:npatch){
    l <- makePatch(context=mtx, spt=pts[np], size=size[np], bgr=bgr, edge=edge, val=val, rast = FALSE)
    if(edge==TRUE){
      eg <- l[[2]]
      l <- l[[1]]
      lst[[np]] <- list(l, eg)
    } else {
      lst[[np]] <- l
    }
    .assignValues(val, l, mtx) #mtx[l] <- val
  }
  if(rast == TRUE) {
    context[unlist(lst)] <- val
    return(context)
  } else {
    return(lst)
  }
}
