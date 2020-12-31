#' Generate a random sibling from a given DNA profile and a given set of allele
#' frequencies
#' 
#' Generates a random sibling from a given DNA profile from a given set of
#' allele frequencies. At each locus, the sibling inherits the first allele of
#' the given profile with one quarter, or the second allele with probability
#' one quarter, both alleles with probability one quarter, or neither with
#' probability one quarter. If the sibling inherits zero or one identical
#' alleles, the missing alleles are chosen at random with probability
#' proportional to the allele frequencies.
#' 
#' The alleles are simply integers rather than the STR repeat numbers. This
#' speeds up computation immensely when calculating any of the LRs or IBS.
#' 
#' @param profile A vector consisting of 2*nLoci elements. Each element in the
#' vector is the (coded) allele held by the individual. This represents the
#' sibling.
#' @param Freqs A list containing two lists labelled loci and freqs. The second
#' list is a list of vectors containing the allele frequencies of each allele
#' at each locus in the multiplex.
#' @return A vector with 2*nLoci elements. Each pair of elements represents the
#' genotpe of the random individual at that locus. The genotype alleles are
#' always ordered so that allele1 <= allele2.
#' @author James M. Curran
#' @seealso randomChild, randomSample
#' @examples
#' 
#' data(fbiCaucs)
#' P1 = randomProfile(fbiCaucs)
#' S1 = randomSib(P1,fbiCaucs)
#' P1
#' S1
#' 
#' @export randomSib
randomSib = function(profile, Freqs){
    nLoci = length(Freqs$loci)
    profSib = rep(0, 2*nLoci)

    for(nLoc in 1:nLoci){
        f = Freqs$freqs[[nLoc]]
        i = sample(1:4, 1)
        a = sample(1:length(f), 2, replace = TRUE, prob = f)
        i1 = 2 * nLoc - 1
        i2 = i1 + 1

        switch(i,
               {profSib[i1:i2] = profile[i1:i2]},
               {profSib[i1:i2] = c(profile[i1], a[1])},
               {profSib[i1:i2] = c(a[1], profile[i2])},
               {profSib[i1:i2] = a}
               )


        if(profSib[i1] > profSib[i2]){
            swap = profSib[i1]
            profSib[i1] = profSib[i2]
            profSib[i2] = swap
        }
    }
    class(profSib) = "profile"
    return(profSib)
}
