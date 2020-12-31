#' Mode
#' 
#' Function calculating the mode.
#' 
#' @param x vector to get mode for
#' @param multimodal wether or not all modes should be returned in case of more than one
#' @param warn should the function warn about multimodal outcomes?
#' 
#' @return vector of mode or modes
#' 
#' @export
#' 
stats_mode <- function(x, multimodal=FALSE, warn=TRUE) {
  res <- stats_mode_multi(x)
  if( identical(multimodal, TRUE) ){
    return(res)
  }else{
    if( warn & length(res) > 1 ){
      warning("modus : multimodal but only one value returned (use warn=FALSE to turn this off)")
    }
    if( !identical(multimodal, FALSE) & length(res) > 1 ){
      return(multimodal)
    }else{
      return(res[1])
    }
  }
}



