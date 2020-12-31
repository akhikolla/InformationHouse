#' Calculate locuswise likelihood ratios for two person victim/suspect mixtures
#' 
#' Calculates the likelihood ratio for pairs of profiles under the propositions
#' \eqn{H_p:\qquad V+S}{Hp: V+S} and \eqn{H_d:\qquad }{Hd: V+U}\eqn{ V+U}{Hd:
#' V+U}, where \eqn{V}{V}, \eqn{S}{S} and \eqn{U}{U} are the victim, the
#' suspect and someone unrelated to the suspect respectively. The calculation
#' does not employ \eqn{\theta}{theta} so there are no assumptions about the
#' subpopulations of the contributors.
#' 
#' 
#' @param profiles A vector of profile lists, from \code{randomProfilePairs}.
#' \code{randomPCPairs} and \code{randomSibPairs} also work but should not
#' really be used as the calculations do not take account of the relationship
#' between the two individuals.
#' @param Freqs A list containing elements \code{freqs}, \code{loci} and
#' \code{counts}. The element \code{freqs} is a list of vectors of allele
#' frequencies at the loci listed in \code{loci}. These frequencies are used to
#' evaluate the LR
#' @return A matrix of LRs calculated at each locus for every pair of profiles.
#' Note this is the set of \eqn{N}{N} profile pairs supplied in
#' \code{profiles}, not a pairwise comparison.
#' @author James M. Curran
#' @examples
#' 
#' data(USCaucs)
#' p = randomProfilePairs(USCaucs, 10000)
#' log.lrs = log10(lrMix(p, USCaucs))
#' boxplot(log.lrs, las = 2)
#' 
#' @export lrMix
lrMix = function(profiles, Freqs){
    N  = length(profiles)
    nLoci = length(Freqs$loci)
    results = matrix(0, nrow = N, ncol = nLoci)

    for(i in 1:N){
      results[i,] = .LRmix(profiles[[i]][[1]],
                        profiles[[i]][[2]],
                        Freqs$freqs)
    }

    colnames(results) = Freqs$loci
    return(results)
}
