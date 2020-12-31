#' Create a linear patch (beta version).
#'
#' @description Create a linear patch, setting direction and convolution. The higher the convolution degree, the weaker the
#' linear shape (and direction).
#' % ADD set length instead of size
#' % FIX values of direction at 0 and 360
#' % FIX values of convol of 1 and 0
#' @inheritParams makePatch
#' @param direction numeric. The direction towards which the patch will point and grow, expressed in degrees between 0 and 360.
#' @param convol numeric. Level of convolution to be assigned to the patch (default is \code{convol=0.5}). A value
#' between 0 and 1, where 0 is no convolution at all (basically a straight line) and 1 is maximum convolution (i.e. tends to form a circular patch).
#' @return A vector of raster cell numbers, or a RasterLayer object if \code{rast=TRUE}. If \code{edge=TRUE} a
#' list of two vectors is returned: one for the inner raster cells and the second for cells at the edge of the patch.
#' @details For any values of \code{convol} > 0.8, no big differences are observed noted. Also direction is progressively lost
#' as convolution increases.
#' @examples
#' library(raster)
#' r <- matrix(0,33,33)
#' r <- raster(r, xmn=0, xmx=10, ymn=0, ymx=10)
#' plot(makeLine(r, size=50, spt = 545, direction=45, convol=0.05, rast=TRUE))
#'
#' @export
makeLine <- function(context, size, direction=NULL, convol=0.5, spt=NULL, bgr=0, edge=FALSE, rast=FALSE, val=1) {
  if(!is.matrix(context)) {
    if(!is.null(spt)){
      if(length(spt) == 2){
        spt <- matrix(spt, ncol=2)
      }
      spt <- .toCellIndex(context, spt)
    }
    mtx <- t(raster::as.matrix(context))
  } else {
    mtx <- context
  }
  bgrCells <- which(mtx == bgr)
  if(is.null(spt)){
    spt <- sample(bgrCells, 1)
  }
  if(length(bgrCells) == 0){
    stop('No background cells available, landscape full. Try checking argument "bgr".')
  }
  if(size > length(bgrCells)){
    warning('Patch size bigger than available landscape.')
  }
  if(.subset(mtx, spt) != bgr){
    wp <- spt
    if(length(bgrCells) > 1){
      spt <- sample(bgrCells, 1)
    } else {
      spt <- bgrCells
    }
    warning('Seed point  ', wp, '  outside background. Re-sampled randomly inside it. New seed:  ', spt)
  }
  .assignValues(val, spt, mtx) # mtx[spt] <- val
  edg <- spt
  cg = 1
  if(!is.null(direction)){
    SD <- .directionDistrib(direction, convol)
  }
  while(cg < size){
    ad <- .contigCells(spt, bgr, mtx)
    if(length(ad) == 0) {
      edg <- edg[edg != spt]
      if(length(edg) <= 1) {
        warning('Maximum patch size reached: ', cg)
        break
      }
      spt <- sample(edg, 1) ## remove to avoid branching (will need more fixes though)
    } else {
      .assignValues(val, ad, mtx) # mtx[ad] <- val
      edg <- c(edg[edg != spt], ad)
      cg <- cg + length(ad)
      if(is.null(direction) | length(ad) == 1){
        spt <- ifelse(length(ad) == 1, ad, sample(ad, 1))
      } else {
        print
        spt <- .direct(spt, ad, direction, SD)
      }
    }
  }
  if(rast == TRUE) {
    context[] <- t(mtx)
    return(context)
  } else if (edge == TRUE) {
    edgVal <- ifelse(val+1 == bgr, val+2, val+1)
    .assignValues(edgVal, edg, mtx) # mtx[edg] <- edgVal
    edg <- which(mtx == edgVal)
    idx <- which(mtx == val)
    return(list(inner = idx, edge = edg))
  } else {
    return( which(mtx == val))
  }
}

#########
.makeString <- function(context, size, direction=NULL, convol=0.5, spt=NULL, bgr=0, edge=FALSE, rast=FALSE, val=1) {
  if(!is.matrix(context)) {
    if(!is.null(spt)){
      if(length(spt) == 2){
        spt <- matrix(spt, ncol=2)
      }
      spt <- .toCellIndex(context, spt)
    }
    mtx <- t(raster::as.matrix(context))
  } else {
    mtx <- context
  }
  bgrCells <- which(mtx == bgr)
  spt <- ifelse(is.null(spt), sample(bgrCells, 1), spt)
  if(length(bgrCells) == 0){
    stop('No background cells available, landscape full. Try checking argument "bgr".')
  }
  if(size > length(bgrCells)){
    warning('Patch size bigger than available landscape.')
  }
  if(.subset(mtx, spt) != bgr){
    wp <- spt
    if(length(bgrCells) > 1){
      spt <- sample(bgrCells, 1)
    } else {
      spt <- bgrCells
    }
    warning('Seed point  ', wp, '  outside background. Re-sampled randomly inside it. New seed:  ', spt)
  }
  edg <- spt
  cg = 1
  if(!is.null(direction)){
    SD <- .directionDistrib(direction, convol)
  }
  while(cg < size){
    ad <- .contigCells(spt, bgr, mtx)
    if(length(ad) == 0) {
      edg <- edg[edg != spt]
      if(length(edg) <= 1) {
        warning('Maximum patch size reached: ', cg)
        break
      }
      spt <- sample(edg, 1) ## remove to avoid branching (will need more fixes though)
    } else {
      .assignValues(val, spt, mtx) # mtx[spt] <- val
      edg <- c(edg[edg != spt], ad)
      cg <- cg + 1
      if(is.null(direction) | length(ad) == 1){
        spt <- ifelse(length(ad) == 1, ad, sample(ad, 1))
      } else {
        print
        spt <- .direct(spt, ad, direction, SD)
      }
    }
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
    .assignValues(edgVal, edg, mtx) # mtx[edg] <- edgVal
    edg <- which(mtx == edgVal)
    idx <- which(mtx == val)
    return(list(inner = idx, edge = edg))
  } else {
    return( which(mtx == val))
  }
}

#########
.directionDistrib <- function(direction, convolution){
  if(direction < 0 | direction > 360){
    stop('direction must be a value between 0 and 360')
  }
  if(convolution < 0 | convolution > 1){
    stop('convolution must be a value between 0 and 1')
  }
  #if(convolution == 1)
  smpSD <- seq(0, 180, 0.2)
  dlt <- sapply(smpSD, function(x){
    q <- stats::qnorm(convolution / 2, direction, x)
    lowQ <- direction - 45
    return(abs(lowQ - q))
  })
  .subset(smpSD, which.min(dlt))
}

.direct <- function(spt, ad, degree, SD){
  ## highly hardcoded on order of adjacent cells in ad:
  ##  # 3 #
  ##  1   2
  ##  # 4 #
  #identify available cells and direction
  dirAd <- sapply(ad, function(x){
    if(x < spt){
      ifelse(x == spt-1, 1, 3)
    } else {
      ifelse(x == spt+1, 2, 4)
    }
  })
  # Sample direction from distribution
  dg <- stats::rnorm(1, degree, SD)
  dg <- ifelse(dg < 0, dg + (trunc(abs(dg/360)) + 1) * 360, dg - 360 * trunc(dg/360))
  dr <- .cellDir(dg, directions=4)
  if(!any(dr %in% dirAd)){  #if sampled direction not available select the closest
    if(dr == 1){
      if(dg >= 270){ dr <- c(3, 4, 2)[c(3, 4, 2) %in% dirAd][1]
      } else { dr <- c(4, 3, 2)[c(4, 3, 2) %in% dirAd][1] }
    } else if(dr == 2){
      if(dg >= 90){ dr <- c(3, 4, 1)[c(3, 4, 1) %in% dirAd][1]
      } else { dr <- c(4, 3, 1)[c(4, 3, 1) %in% dirAd][1] }
    } else if(dr == 3){
      if(dg < 45){ dr <- c(2, 1, 4)[c(2, 1, 4) %in% dirAd][1]
      } else { dr <- c(1, 2, 4)[c(1, 2, 4) %in% dirAd][1] }
    } else {
      if(dg >= 180){ dr <- c(1, 2, 3)[c(1, 2, 3) %in% dirAd][1]
      } else { dr <- c(2, 1, 3)[c(2, 1, 3) %in% dirAd][1] }
    }
  }
  ad[dirAd == dr]
}

.cellDir <- function(directed, directions){
  if(directions == 4){
    cellDir <- sapply(directed, function(x){ ## Faster and cleaner than 'cut' function
      if(x < 45 | x >= 315){
        3
      } else if(x >= 45 & x < 135){
        2
      } else if(x >= 135 & x < 225){
        4
      } else {
        1
      }
    })
  } else {
    cellDir <- sapply(directed, function(x){
      if(x >= 337.5 | x < 22.5){
        1
      } else if(x >= 22.5 & x < 67.5){
        2
      } else if(x >= 67.5 & x < 112.5){
        3
      } else if(x >= 112.5 & x < 157.5){
        4
      } else if(x >= 157.5 & x < 202.5){
        5
      } else if(x >= 202.5 & x < 247.5){
        6
      } else if(x >= 247.5 & x < 292.5){
        7
      } else {
        8
      }
    })
  }
  cellDir
}

#########
.direction <- function(degree, convolution=1, size, directions=4){
  if(degree < 0 | degree > 360){stop('degree must be a value between 0 and 360')}
  if(convolution < 0 | convolution > 1){stop('convolution must be a value between 0 and 1')}
  angle <- ifelse(directions==4, 90, 45) / 2
  smpSD <- seq(0, 180, 0.2)
  dlt <- sapply(smpSD, function(x){
    q <- stats::qnorm(convolution/2, degree, x)
    lowQ <- degree - angle
    return(abs(lowQ - q))
  })
  SD <- smpSD[which.min(dlt)]
  directed <- stats::rnorm(size-1, degree, SD)
  directed <- ifelse(directed < 0, directed + (trunc(abs(directed/360)) + 1) * 360, directed - 360 * trunc(directed/360))
  .cellDir(directed, directions=directions)
}




