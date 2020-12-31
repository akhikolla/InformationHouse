#' Generate one or more random parent/child pairs from a given set of allele
#' frequencies
#' 
#' Generates one or more pairs random parent/child pairs from a given set of
#' allele frequencies.
#' 
#' The alleles are simply integers rather than the STR repeat numbers. This
#' speeds up computation immensely when calculating any of the LRs or IBS.
#' 
#' @param Freqs A list containing two lists labelled loci and freqs. The second
#' list is a list of vectors containing the allele frequencies of each allele
#' at each locus in the multiplex.
#' @param BlockSize The number of pairs of profiles to generate
#' @return A list of length \code{BlockSize}. Each element of the list has a
#' sublist containing two profiles called \code{parent} and \code{child}
#' @author James M. Curran
#' @seealso randomSibPairs, randomProfilePairs
#' @examples
#' 
#' data(fbiCaucs)
#' P = randomPCPairs(fbiCaucs)
#' P$parent
#' P$child
#' 
#' @export randomPCPairs
randomPCPairs = function(Freqs, BlockSize = 1){
    nLoci = length(Freqs$loci)
    Profile = vector(mode = "list", length = BlockSize)

    Parent = .randomProfiles(Freqs$freqs, BlockSize)
    Child = randomChildren(Parent, Freqs$freqs, BlockSize)

    for(b in 1:BlockSize){
        i1 = (b - 1) * 2 * nLoci + 1
        i2 =  b * 2 * nLoci

        Profile[[b]] = vector(mode = "list", length = 2)
        names(Profile[[b]]) = c("parent", "child")

        Profile[[b]]$parent = Parent[i1:i2]
        class(Profile[[b]]$parent) = "profile"

        Profile[[b]]$child = Child[i1:i2]
        class(Profile[[b]]$child) = "profile"
    }

    if(BlockSize == 1){
        return(Profile[[1]])
    }else{
        return(Profile)
    }
}
