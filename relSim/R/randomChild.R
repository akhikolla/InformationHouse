#' Generate a random child from a given DNA profile and a given set of allele
#' frequencies
#' 
#' Generates a random child (or parent) from a given DNA profile from a given
#' set of allele frequencies. At each locus, the child inherits the first
#' allele of the given profile with one half, or the second allele with
#' probability one half. The second allele is chosen at random with probability
#' proportional to the allele frequencies.
#' 
#' The alleles are simply integers rather than the STR repeat numbers. This
#' speeds up computation immensely when calculating any of the LRs or IBS.
#' 
#' @param profile A vector of length 2*nLoci. Each entry in the vector is the
#' (coded) allele held by the individual. This represents the parent. The
#' relationship is reflexive so it does not matter if the profile is a parent
#' or a child.
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
#' C1 = randomChild(P1,fbiCaucs)
#' P1
#' C1
#' 
#' @importFrom stats runif
#' @export randomChild
randomChild = function(profile, Freqs){
    nLoci = length(Freqs$loci)
    profChild = rep(0, 2*nLoci)

    for(nLoc in 1:nLoci){
        f = Freqs$freqs[[nLoc]]
        a = sample(1:length(f), 1, prob = f)
        u = runif(1)
        i1 = 2*nLoc - 1
        i2 = i1 + 1

        if(u < 0.5){
            profChild[i1:i2] = c(profile[i1], a)
        }else{
            profChild[i1:i2] = c(a, profile[i2])
        }

        if(profChild[i1] > profChild[i2]){
            swap = profChild[i1]
            profChild[i1] = profChild[i2]
            profChild[i2] = swap
        }
    }

    class(profChild) = "profile"

    return(profChild)
}
