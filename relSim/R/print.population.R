#' Print summary details of a substructed population
#' 
#' Nicely prints summary information about a substructured population created
#' using \code{breedFst}
#' 
#' 
#' @param x The population object to be printed
#' @param \dots Ignored - really should be passed to print, but given cat is
#' actually called they are ignored
#' @author James M. Curran
#' @seealso breedFst
#' @examples
#' 
#' data(fbiCaucs)
#' p = breedFst(fbiCaucs)
#' print(p)
#' 
print.population = function(x,...){
    cat("Substructured population\n")
    cat("------------------------\n")
    cat(paste("Number of members:", x$nProfiles, "\n"))
    cat(paste("Number of subpopulations", x$nSubpops, "\n"))
    cat(paste("Number of loci:", x$nLoci, "\n"))
    cat(paste("Expected level of coancestry (theta):", x$theta, "\n"))
}
