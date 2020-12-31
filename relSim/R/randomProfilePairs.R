#' Generate one or more random DNA profile pairs from a given set of allele
#' frequencies
#' 
#' Generates one or more random DNA profile pairs from a given set of allele
#' frequencies.
#' 
#' The alleles are simply integers rather than the STR repeat numbers. This
#' speeds up computation immensely when calculating any of the LRs or IBS.
#' 
#' @param Freqs A list containing two lists labelled loci and freqs. The second
#' list is a list of vectors containing the allele frequencies of each allele
#' at each locus in the multiplex.
#' @param BlockSize The number of pairs of profiles to generate
#' @return A list of length \code{BlockSize}. Each element of the list has a
#' sublist containing two profiles called \code{prof1} and \code{prof2}
#' @author James M. Curran
#' @seealso randomPCPairs, randomSibPairs
#' @examples
#' 
#' data(fbiCaucs)
#' P = randomProfilePairs(fbiCaucs)
#' P$prof1
#' P$prof2
#' 
#' @export randomProfilePairs
randomProfilePairs = function(Freqs, BlockSize = 1){
    ## nLoci = length(Freqs$loci)
    ## profile = matrix(0, nc = 2, nr = nLoci)

    ## for(nLoc in 1:nLoci){
    ##     f = Freqs$freqs[[nLoc]]
    ##     profile[nLoc,] = sample(1:length(f), 2, replace = T, prob = f)

    ##     if(profile[nLoc,1] > profile[nLoc, 2]){
    ##         swap = profile[nLoc, 1]
    ##         profile[nLoc, 1] = profile[nLoc, 2]
    ##         profile[nLoc, 2] = swap
    ##     }
    ## }

    ## class(profile) = "profile"
    ## return(profile)

    

    prof1 = .randomProfiles(Freqs$freqs, BlockSize)
    prof2 = .randomProfiles(Freqs$freqs, BlockSize)
    nLoci = length(Freqs$loci)
    Profile = vector(mode = "list", length = BlockSize)
    
    for(b in 1:BlockSize){
        i1 = (b - 1) * 2 * nLoci + 1
        i2 =  b * 2 * nLoci

        Profile[[b]]$prof1 = prof1[i1:i2]
        Profile[[b]]$prof2 = prof2[i1:i2]
        class(Profile[[b]]$prof1) = "profile"
        class(Profile[[b]]$prof2) = "profile"
    }

    if(BlockSize==1){
        return(Profile[[1]])
    }else{
        return(Profile)
    }
}
