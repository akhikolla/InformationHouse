#' Generate a random DNA profile from a given set of allele frequencies
#' 
#' Generates a random DNA profile from a given set of allele frequencies.
#' 
#' The alleles are simply integers rather than the STR repeat numbers. This
#' speeds up computation immensely when calculating any of the LRs or IBS.
#' 
#' @param Freqs A list containing two lists labelled loci and freqs. The second
#' list is a list of vectors containing the allele frequencies of each allele
#' at each locus in the multiplex.
#' @return A vector with 2*nLoci elements. Each pair of elements represents the
#' genotpe of the random individual at that locus. The genotype alleles are
#' always ordered so that allele1 <= allele2.
#' @author James M. Curran
#' @seealso randomChild, randomSample, randomSib
#' @examples
#' 
#' data(fbiCaucs)
#' P1 = randomProfile(fbiCaucs)
#' 
#' @export randomProfile
randomProfile = function(Freqs){
    nLoci = length(Freqs$loci)
    profile = rep(0, 2*nLoci)

    for(nLoc in 1:nLoci){
        i1 = 2*nLoc - 1
        i2 = i1 + 1
        f = Freqs$freqs[[nLoc]]
        profile[i1:i2] = sample(1:length(f), 2, replace = T, prob = f)

        if(profile[i1] > profile[i2]){
            swap = profile[i1]
            profile[i1] = profile[i2]
            profile[i2] = swap
        }
    }

    class(profile) = "profile"
    return(profile)

}
