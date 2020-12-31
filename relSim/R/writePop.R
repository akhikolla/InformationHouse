#' Saves/writes population profiles to disk
#' 
#' Writes a population of profiles to disk using the original allele
#' designations rather than the internal integer representations that are used
#' for the other functions.
#' 
#' 
#' @param fileName the name and path where the population profiles are to be
#' saved to.
#' @param pop an object of class \code{population}, most likely produced by
#' \code{breedFst}
#' @param addAmelo The simulated populations do not have Amelogenin. If
#' \code{TRUE} then an Amelogenin marker is added to the population, and all
#' the profiles are set to male XY, although this is coded to 1,2 to keep the
#' allele designations numeric.
#' @param delim The allele delimiter.
#' @param dupLoci If \code{TRUE} the locus names are written twice in the header, otherwise just once.
#' @note Rare alleles are recoded to 108.1.
#' @author James M. Curran
#' @seealso breedFst
#' @examples
#' 
#' data(USCaucs)
#' pop = breedFst(USCaucs)
#' \dontrun{
#'   writePop("USCaucs.csv", pop)
#'   }
#' @export
writePop = function(fileName,  pop, addAmelo = FALSE, delim = ',', dupLoci = TRUE){
  nLoci = pop$nLoci
  Alleles = lapply(pop$Freqs$freqs, names)
  Alleles = lapply(Alleles, function(x){x[grep("R", x)] = "108.1"; x})
  
  toProf = function(prof){
    A = rep("", 2 * nLoci)
    if(addAmelo){
      A = c(A, "1", "2")
    }
    
    for(i in 1:nLoci){
      i1 = 2 * i - 1
      i2 = i1 + 1
      A[i1] = Alleles[[i]][prof[i1]]
      A[i2] = Alleles[[i]][prof[i2]]
    }
    
    return(paste(A, collapse = delim))
  }
  
  popMatrix = matrix(pop$profiles, nrow = pop$nProfiles, byrow = T) 
  profiles = apply(popMatrix, 1, toProf)
  f1 = file(fileName, 'w')
  locusHeader = rep(pop$Freqs$loci, rep(2, nLoci))
  if(addAmelo)
    locusHeader = paste0(locusHeader, ',"Amelo", "Amelo"')
  writeLines(paste0('"', paste(rep(pop$Freqs$loci, rep(2, nLoci)), collapse = '","'), '"'), f1)
  writeLines(profiles, f1)
  close(f1)
}
