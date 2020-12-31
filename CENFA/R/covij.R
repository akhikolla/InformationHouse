#' @keywords internal
# @importFrom raster cellStats filename extension

# This function calculates the covariance between two raster layers.

.covij <- function(x, y, w, sample = sample){

  sm <- canProcessInMemory(x)

  if(sm){
    x <- values(x)
    x <- x - mean(x, na.rm = T)
    y <- values(y)
    y <- y - mean(y, na.rm = T)
    if(!is.null(w)){
      w <- values(w)
      sumw <- sum(w, na.rm = T)
      #if(sumw > 1) w <- w / (sumw)
      r <- na.omit(w * x * y)
      v <- sum(r, na.rm = T)/(sumw - sample)
    } else {
      r <- na.omit(x * y)
      nn <- length(r)
      v <- sum(r, na.rm = T)/(nn - sample)
    }
  } else if(!sm){

    x <- scale(x, scale = F)
    y <- scale(y, scale = F)
    if(!is.null(w)){
      sumw <- cellStats(w, sum, na.rm = T)
      #if(sumw > 1) w <- w / (sumw)
      r <- w * x * y
      v <- cellStats(r, stat = 'sum', na.rm = T) / (sumw - sample)
    } else {
      r <- x * y
      nn <- length(r[!is.na(r)])
      v <- cellStats(r, stat = 'sum', na.rm = T) / (nn - sample)
    }
    f <- filename(r)
    file.remove(c(f, extension(f, '.gri')))
  }
  return(v)
}
