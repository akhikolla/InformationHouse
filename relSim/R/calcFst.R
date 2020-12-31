#' Caculate locus-wise and population \eqn{F_{ST}}{Fst} values
#' 
#' This procedure uses the method of Weir and Cockerham to estimate
#' \eqn{\theta}{theta} (\eqn{F_{ST}}{Fst}) for a population with substructure
#' 
#' 
#' @param Pop An object type 'population'
#' @param subPopIdx If this vector is not null, then it must consist of
#' \eqn{N}{N} elements with values from 1 to \eqn{n_s}{ns} representing which
#' subpopulation each member of \code{Pop$profiles} belongs to. If it is null
#' then it is assumed that the population consists of \eqn{n_s}{ns}
#' subpopulations of equal size \eqn{N_s}{Ns} so that \eqn{n_s\times N_s =
#' N}{ns*Ns = N}
#' @return A vector of length \eqn{n_{loci}+1}{nloci+1} with locus-wise
#' \eqn{\theta}{theta} values and an overall \eqn{\theta}{theta} value for the
#' population
#' @author James M. Curran
#' @seealso breedFst
#' @references Weir, B.S., Genetic Data Analysis II, (1996) p.173--179,
#' Sinauer, Sunderland, MA.
#' @examples
#' 
#' data(USCaucs)
#' p = breedFst(USCaucs)
#' fst = calcFst(p)
#' fst
#' 
#' @export calcFst
calcFst = function(Pop, subPopIdx = NULL){
    if(class(Pop) != "population")
        stop("Pop must be of class Population")

    nLoci = Pop$nLoci
    NumLocusAlleles = sapply(Pop$Freqs$freqs, length)

    if(is.data.frame(Pop$profiles)){
        Pop$profiles = as.vector(t(as.matrix(Pop$profiles[,-1])))
    }

    N = length(Pop$profiles)

    if(N %% nLoci !=0){
        stop("This proceedure only works for complete profiles")
    }

    N = Pop$nProfiles
    ns = Pop$nSubpops

    if(N %% ns != 0)
        stop("The number of subpopulations must evenly divide N")

    if(is.null(subPopIdx)){
        Ns = N / ns
        subPopIdx = rep(1:ns, rep(Ns, ns))
    }else{
        if(length(subPopIdx != N))
            stop("subPopIdx length must equal pop length")

        ns = length(unique(subPopIdx))
    }

    Fst = .calcFst(Pop = Pop$profiles,
                   SubPopIdx = subPopIdx,
                   N = N, 
                   ns = ns,
                   nLoci = nLoci,
                   NumLocusAlleles = NumLocusAlleles)

    names(Fst)[1:nLoci] = Pop$Freqs$loci
    names(Fst)[nLoci + 1] = "Overall"
    return(Fst)
}
