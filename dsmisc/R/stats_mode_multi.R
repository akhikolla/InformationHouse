#' Mode Allowing for Multi Modal Mode
#'
#' Function calculating the mode, allowing for multiple modes in case of equal frequencies. 
#' 
#'
#' @inheritParams stats_mode
#'
#' @return vector with all modes
#' 
#' @export
#'
stats_mode_multi <- 
  function(x){
    x_unique <- unique(x)
    tab_x    <- tabulate(match(x, x_unique))
    x_unique[which(tab_x==max(tab_x))]
  }
