#' Print a DNA profile
#' 
#' Nicely prints a profile object out in genotype pairs
#' 
#' 
#' @param x The profile object to be printed
#' @param horizontal if \code{TRUE} then the profile will print on a single
#' line instead of multiple lines. Useful for comparing two profiles
#' @param \dots Ignored - really should be passed to print, but given cat is
#' actually called they are ignored
#' @author James M. Curran
#' @examples
#' 
#' data(fbiCaucs)
#' P1 = randomProfile(fbiCaucs)
#' P2 = randomProfile(fbiCaucs)
#' P1
#' print(P1, horizontal = TRUE)
#' print(P2, horizontal = TRUE)
#' 
print.profile = function(x, horizontal = FALSE, ...){
    nLoci = length(x)/2
    i1 = 2*(1:nLoci) - 1
    i2 = i1 + 1
    strProf = sprintf("%2d/%-2d", x[i1], x[i2])

    if(horizontal){
        cat(paste(paste(strProf, collapse = ""), "\n"))
    }else{
        cat(paste(paste(strProf, collapse = "\n"), "\n"))
    }
}
