#' Calculate the exclusion power of a multiplex by locus
#' 
#' Calculates the exclusion power
#' \deqn{1-2\left(\sum_{i=1}^{n_l}p_i^2\right)^2-4\sum_{i=1}^{n_l}p_i^4}{1-2*sum(pi^2)^2-4*sum(pi^4)}
#' at each locus for a set of allele frequencies.
#' 
#' 
#' @aliases exclusionPower ep
#' @param Freqs A list containing two vectors and a list, called loci, counts,
#' and freqs. The elements of loci are the loci present in the multiplex. The
#' elements are freqs a vectors of allele frequencies for the locus. The
#' elements of counts are irrelevant here.
#' @return The exclusion power for each locus.
#' @author James M. Curran
#' @references NRC II, Evaluation of Forensic Evidence, (1996), p.96, National
#' Academy Press.
#' @examples
#' 
#' data(USCaucs)
#' ep(USCaucs)
#' 
#' ## get the multiplex wide exclusion power
#' 1 - prod(1-ep(USCaucs))
#' 
#' @export exclusionPower
exclusionPower = function(Freqs){
    exclPwr = 1 - sapply(Freqs$freqs,function(p){2 * sum(p^2)^2 - sum(p^4)})
    names(exclPwr) = Freqs$loci
    return(exclPwr)
}

#' @export ep
ep = exclusionPower

