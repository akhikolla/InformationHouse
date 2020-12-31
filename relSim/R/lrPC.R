#' Likelihood Ratio for Parent-Child / Paternity Index
#' 
#' Calculates Likelihood Ratio comparing the probability of two profiles if
#' they are indeed parent-child compared to unrelated. This is the paternity
#' index or PI.
#' 
#' 
#' @param parent A matrix consisting of 2 columns and nLoci rows. Each entry in
#' the matrix is the (coded) allele held by the individual. This represents the
#' alleged parent. The relationship is reflexive so it does not matter which
#' profile is labelled parent and child.
#' @param child See \code{parent}
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
#' profiles are parent and child. There is no mutation built into this
#' calculation. This means that the LR will be zero if the profiles do not
#' share at least one allele in common at each locus in the multiplex.
#' @author James M. Curran
#' @seealso lrSib, IBS
#' @references Buckleton, J, Triggs, C.M., and Walsh, S.J. (2005)\emph{Forensic
#' DNA Evidence Interpretation}, CRC Press., Boca Raton, FL. p.410
#' @examples
#' 
#' data(fbiCaucs)
#' P1 = randomProfile(fbiCaucs)
#' C1 = randomChild(P1, fbiCaucs)
#' lrPC(P1, C1, fbiCaucs)
#' 
#' @export lrPC
lrPC = function(parent, child, Freqs = NULL,
                nLoci = length(parent)/2, f = NULL, n = NULL){
  
  if(is.null(Freqs) & (is.null(f) | is.null(n))){
    stop("Must specify either Freqs or both f and n")
  }
  
  if(is.null(Freqs)){
    loc = rep(1:length(n), n)
    Freqs = split(f, n)
  }
  
  
  lr = .lrPC(parent, child, Freqs$freqs)

  return(lr)
}
