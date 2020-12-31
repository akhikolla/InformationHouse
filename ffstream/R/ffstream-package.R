#' ffstream package 0.1.6
#' 
#' This function gives information on the package and how to see the vignette.
#'
ffstream <- function(){
    cat("See the ffstream documentation or vignette for more details.\n")
    cat("To see the command to run the vignette, ")
    cat("run the function `ffstream_vignette()`\n")
}


#' ffstream vignette 0.1.6
#' 
#' This function opens the ffstream vignette.
#'
#' @examples
#' ffstream_vignette()
#' @export
ffstream_vignette <- function(){
    #vignette("intro", package="ffstream") 
    cat('To see the vignette, run `vignette("intro", package="ffstream")`\n')
}
