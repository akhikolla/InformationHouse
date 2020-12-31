#' Likelihood Ratio / Kinship Index for full-siblings
#' 
#' Calculates Likelihood Ratio comparing the probability of two profiles if
#' they are indeed full-sibs compared to unrelated. This is sometimes called
#' the kinship index (KI) for full-sibs.
#' 
#' 
#' @param sib1 A matrix consisting of 2 columns and nLoci rows. Each entry in
#' the matrix is the (coded) allele held by the individual. This represents the
#' alleged sibling. The relationship is reflexive so it does not matter which
#' profile is labelled sib1 and sib2.
#' @param sib2 See \code{sib1}
#' @param Freqs A list containing two lists labelled loci and freqs. The second
#' list is a list of vectors containing the allele frequencies of each allele
#' at each locus in the multiplex. This argument or both f and n must be
#' specified
#' @param nLoci The number of loci in the profiles
#' @param f A concatenated vector of allele frequencies. Specifying this speeds
#' up computation enormously
#' @param n A vector of length \code{nLoci} giving the number of alleles at
#' each locus. Specifying this in advance enormously speeds up computation
#' @return A value between 0 and infinity representing support (or lack of
#' support if the value is less than 1) for the hypothesis that the two
#' profiles are full-siblings. There is no mutation built into this
#' calculation.
#' @author James M. Curran
#' @seealso lrSibDebug, lrPC, IBS
#' @references Buckleton, J, Triggs, C.M., and Walsh, S.J. (2005)\emph{Forensic
#' DNA Evidence Interpretation}, CRC Press., Boca Raton, FL. p.411
#' @examples
#' 
#' data(fbiCaucs)
#' P1 = randomProfile(fbiCaucs)
#' S1 = randomSib(P1, fbiCaucs)
#' P2 = randomProfile(fbiCaucs)
#' lrSib(P1, S1, fbiCaucs)
#' lrSib(P1, P2, fbiCaucs)
#' 
#' @export lrSib
lrSib = function(sib1, sib2, Freqs  = NULL, nLoci = length(sib1)/2,
                 f = NULL, n = NULL){
  
  if(is.null(Freqs) & (is.null(f) | is.null(n))){
    stop("Must specify either Freqs or both f and n")
  }
  
  if(is.null(Freqs)){
    loc = rep(1:length(n), n)
    Freqs = split(f, n)
  }
  
  
  lr = .lrSib(sib1, sib2, Freqs$freqs)
  
  return (lr)
}
