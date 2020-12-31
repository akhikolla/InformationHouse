#' Identity by state
#' 
#' Calculates the total number of alleles that are shared by two profiles. If
#' the two profiles in question are indeed relatives then the matching alleles
#' may be identical by descent, or by random chance alone, hence identity by
#' state.
#' 
#' 
#' @param prof1 A matrix consisting of 2 columns and nLoci rows. Each entry in
#' the matrix is the (coded) allele held by the individual.
#' @param prof2 See \code{prof1}
#' @param nLoci The number of loci in the profiles. Specifying this value
#' speeds up computation enormously.
#' @param bPrint If true then the result is printed locus by locus. This
#' feature exists primarily for debugging purposes.
#' @return An integer between 0 and 2*nLoci representing the total number of
#' alleles that match in the two profiles.
#' @author James M. Curran
#' @examples
#' 
#' data(fbiCaucs)
#' P1 = randomProfile(fbiCaucs)
#' C1 = randomChild(P1, fbiCaucs)
#' IBS(P1, C1)
#' IBS(P1, C1, bPrint = TRUE)
#' 
#' @export IBS
IBS = function(prof1, prof2, nLoci = length(prof1) / 2, bPrint = FALSE){
    results =  .IBS_Caller(prof1, prof2, nLoci)

    if(bPrint){
        for(i in 1:nLoci){
            x = matrix(c(prof1[c(2 *i - 1, 2 * i)], prof2[c(2 * i - 1, 2 * i)]), nrow = 1)
            m = locusIBS(x)

            if(bPrint){
                cat(paste(prof1[c(2 * i - 1, 2 * i)], collapse = "/"),
                    "\t",
                    paste(prof2[c(2 * i - 1, 2 * i)], collapse = "/"),
                    "\t",
                    ifelse(m > 0, "TRUE", "FALSE"),
                    "\t",
                    m,"\n", sep = "")
            }
        }
    }

    return(results)
}
