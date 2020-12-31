#' Create a single patch
#'
#' @description Function will create a single patch. \strong{NOTE}: \code{makeClass} should be used preferably when creating a single patch, as better error and exception handling is provided.
#' @param context Raster object or matrix, an empty landscape raster or a mask indicating where the patch cannot be generated (see bgr below).
#' @param size integer. Size of the patch to be generated, as number of raster cells.
#' @param spt integer or matrix. The seed point location around which the patch is generated (a random point is given by default). It can be an integer, as index of the cell in the raster, or a two columns matrix indicating x and y coordinates (an integer vector of length 2 is accepted too).
#' @param bgr integer. A single value of background cells, where a patch can be generated (default is zero). Cells/classes which cannot be changed must have a different value.
#' @param edge logical. Should the vector of edge cells of the patch be returned?
#' @param rast logical. If TRUE returns a Raster object, otherwise a vector of cell numbers where the patch occurs
#' @param val integer. The value to be assigned to patch cells, when \code{rast=TRUE}
#' @return A vector of raster cell numbers, or a RasterLayer object if \code{rast=TRUE}. If \code{edge=TRUE} a
#' list of two vectors is returned: one for the inner raster cells and the second for cells at the edge of the patch.
#' @details The patch is created starting from the seed point and iteratively sampling randomly neighbouring cells at the edge of the patch, according to von Neumann neighbourhood (four cells, aka Rook case).
#' There is a tolerance of +/- 3 cells from the patch size declared in \code{size} argument.
#' Argument \code{bgr} accepts a single value only, unlike \code{makeClass} that accepts multiple and should therefore preferred.
#' @examples
#' library(raster)
#' mtx = matrix(0, 33, 33)
#' r = raster(mtx, xmn=0, xmx=10, ymn=0, ymx=10)
#' patchSize = 500
#' rr = makePatch(r, patchSize, rast=TRUE)
#' plot(rr)
#'
#' ## Create a patch with value 3, starting from the centre cell
#' mtx = matrix(0, 33, 33)
#' r = raster(mtx, xmn=0, xmx=10, ymn=0, ymx=10)
#' newVal = 3
#' centre = 545
#' cells = makePatch(r, patchSize, centre)
#' r[cells] = newVal
#' plot(r)
#'
#' ## Now create a new patch with value 10 and size 100 inside the existing patch
#' rr = makePatch(r, 100, bgr=newVal, rast=TRUE, val=10)
#' plot(rr)
#' @export
makePatch <- function(context, size, spt=NULL, bgr=0, edge=FALSE, rast=FALSE, val=1) {
  if(!is.matrix(context)) {
    if(!is.null(spt)){
      if(length(spt) == 2){
        spt <- matrix(spt, ncol=2)
      }
      spt <- .toCellIndex(context, spt)
    }
    mtx <- t(raster::as.matrix(context))
    warningSwitch <- TRUE
  } else {
    mtx <- context
    warningSwitch <- FALSE
  }
  bgrCells <- which(mtx == bgr)
  if(length(bgrCells) == 0){
    stop('No background cells available with value ', bgr, '. Try checking argument "bgr".')
  }
  if(size > length(bgrCells)){
    warning('Patch size bigger than available background.')
  }
  if(is.null(spt)){
    spt <- sample(bgrCells, 1)
  }
  if(spt > length(context) | spt < 1 | spt %% 1 != 0){
    stop('Seed point not valid. Must be an integer between 1 and the total number of cells of "context".')
  }
  if(.subset(mtx, spt) != bgr){ #mtx[spt] != bgr
    wp <- spt
    if(length(bgrCells) > 1){
      spt <- sample(bgrCells, 1)# bgrCells[.Internal(sample(length(bgrCells), 1L, FALSE, NULL))] #
    } else {
      spt <- bgrCells
    }
    if(warningSwitch){
      warning('Seed point  ', wp, '  outside background. Re-sampled randomly inside it. New seed:  ', spt)
    }
  }
  .assignValues(val, spt, mtx) #mtx[spt] <- val
  edg <- spt
  cg = 1
  while(cg < size){
    ad <- .contigCells(spt, bgr, mtx)
    if(length(ad) == 0) {
      edg <- edg[edg != spt]
      if(length(edg) <= 1) {
        warning('Patch size reached from seed point ', spt, ' was ', cg, ' . No further background cells available for the patch.')
        break
      }
      spt <- sample(edg, 1) # edg[.Internal(sample(length(edg), 1L, FALSE, NULL))] #
    } else {
      .assignValues(val, ad, mtx) #mtx[ad] <- val
      edg <- c(edg[edg != spt], ad)
      cg <- cg + length(ad)
      if(length(edg) == 1){
        spt <- edg
      } else {
        spt <- sample(edg, 1) # edg[.Internal(sample(length(edg), 1L, FALSE, NULL))] #
      }
    }
    #Rcpp::checkUserInterrupt()
  }
  if(rast == TRUE) {
    context[] <- t(mtx)
    return(context)
  } else if (edge == TRUE) {
    if(val+1 == bgr){
      edgVal <- val+2
    } else {
      edgVal <- val+1
    }
    .assignValues(edgVal, edg, mtx) #mtx[edg] <- edgVal
    edg <- which(mtx == edgVal)
    idx <- which(mtx == val)
    return(list(inner = idx, edge = edg))
  } else {
    return( bgrCells[.subset(mtx, bgrCells) == val] ) #bgrCells[mtx[bgrCells] == val]
  }
}

## Converts matrix of coordinates into cell number
.toCellIndex <- function(rast, coord) {
  if(is.vector(coord)){
    coord
  } else if(is.matrix(rast)){
    stop('Cannot use spatial coordinates on a matrix. Please provide a raster or a vector of cell indexes.')
  } else {
    raster::cellFromXY(rast, coord)
  }
}

