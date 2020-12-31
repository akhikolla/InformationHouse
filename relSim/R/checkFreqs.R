#' Make sure that the frequencies are such
#' 
#' Checks whether a list of frequencies at a series of genetic loci both sum to
#' one and lie between 0 and 1.
#' 
#' If a locus fails to sum to one, or there are alleles which fall below zero
#' or above one, then a warning message will be returned for each item in
#' error.
#' 
#' @param Freqs A list containg elements \code{loci} and \code{freqs}.
#' \code{freqs} is a list of vectors containing the frequencies at the given
#' loci.
#' @author James M. Curran
#' @seealso normalizeFreqs
#' @examples
#' 
#' data(fbiCaucs)
#' checkFreqs(fbiCaucs)
#' 
#' ## induce an error
#' fbiCaucs$freqs[[1]] = runif(10)
#' checkFreqs(fbiCaucs)
#' 
#' @export checkFreqs
checkFreqs = function(Freqs){
    sums = sapply(Freqs$freqs, sum)
    ranges = sapply(Freqs$freqs, function(x){any(x < 0 | x > 1)})

    if(any(sums!=1)){
        Loci = Freqs$loci[sums!=1]
        sums = sums[sums!=1]
        cat(paste("The locus", paste(Loci,sums), "does not sum to 1\n"))
    }

    if(any(ranges)){
        Loci = Freqs$loci[ranges]
        cat(paste("The loci", paste(Loci),
                  "have values less than 0 or greater than 1\n"))
    }
}
