#' Identity by state at a locus
#' 
#' Calculates the number of alleles that are shared by two profiles at a single
#' locus.  If the two profiles in question are indeed relatives then the
#' matching alleles may be identical by descent, or by random chance alone,
#' hence identity by state.
#' 
#' 
#' @param profMat A matrix consisting of 4 columns and N rows. Each row in the
#' matrix consists of the genotypes of two individuals.
#' @return A vector of length N containing values 0, 1, or 2 depending on how
#' many alleles each pair of profiles share at a locus.
#' @author James M. Curran
#' @examples
#' 
#' data(fbiCaucs)
#' G = randomSample(1, fbiCaucs, rel = 'FS', N = 1000)
#' ibs = locusIBS(G)
#' barplot(tabulate(ibs+1, nbins = 3))
#' 
#' @export locusIBS
locusIBS = function(profMat){

    ## profIbs = function(prof1, prof2){
    ##     res = 0

    ##     a1 = prof1[1]
    ##     a2 = prof1[2]
    ##     b1 = prof2[1]
    ##     b2 = prof2[2]

    ##     if(a1 == b1 && a2 == b2){
    ##         res = 2
    ##     }else if((a1 == b1) || (a2 == b1) || (a1 == b2) || (a2 == b2)){
    ##         res = 1
    ##     }else{
    ##         res = 0
    ##     }

    ##     return(res)
    ## }

    N = nrow(profMat)
    nc = ncol(profMat)

    if(nc != 4)
        stop("Wrong dimensions")

    p = as.vector(t(profMat))
    return(.locusIBS(p, N))
}
