#' Normalize frequencies to 1
#' 
#' Normalize a list of frequencies at a series of genetic loci both sum to one.
#' Not that this does not deal with the problem of values larger than one or
#' smaller than zero.
#' 
#' Divides vector in Freqs$freqs by the vector sum.
#' 
#' @param Freqs A list containg elements \code{loci} and \code{freqs}.
#' \code{freqs} is a list of vectors containing the frequencies at the given
#' loci.
#' @return A list containg elements \code{loci} and \code{freqs}. \code{freqs}
#' is a list of vectors containing the frequencies at the given loci.
#' @author James M. Curran
#' @seealso checkFreqs
#' @examples
#' 
#' data(fbiCaucs)
#' 
#' ## induce an error
#' fbiCaucs$freqs[[1]] = rgamma(10,1,1)
#' checkFreqs(fbiCaucs)
#' fbiCaucs = normalizeFreqs(fbiCaucs)
#' checkFreqs(fbiCaucs)
#' 
#' @export normalizeFreqs
normalizeFreqs = function(Freqs){
    nLoci = length(Freqs$loci)

    for(i in 1:nLoci){
        Freqs$freqs[[i]] = Freqs$freqs[[i]] / sum(Freqs$freqs[[i]])
    }

    return(Freqs)
}
