#' Saves/writes population frequencies to disk
#' 
#' Writes a population of profiles to disk using the original allele
#' designations rather than the internal integer representations that are used
#' for the other functions.
#' 
#' 
#' @param fileName the name and path where the population profiles are to be
#' saved to.
#' @param pop a \code{list} containing elements \code{loci} and \code{freqs}. \code{loci}
#' is a vector with the loci in the data set. \code{freqs} is a list of vectors with elements
#' named after the elements in \code{loci}. Each locus in \code{freqs} is a vector of allele frequencies
#' with the allele names given by the named elements.
#' \code{TRUE} then an Amelogenin marker is added to the population, and all
#' the profiles are set to male XY, although this is coded to 1,2 to keep the
#' allele designations numeric.
#' @param n the number of people in the database. This is arbitrarily set to 100 by default.
#' @param delim The allele delimiter.
#' @return a matrix which contains the table written to file.
#' @note Rare alleles are recoded to 108.1. This is unlikely to do the right thing
#' when you have things like <5 or >20 in your allele names. Given it is impossible to 
#' predict what a user would like to do, I suggest you recode them yourself before using
#' this function.
#' @author James M. Curran
#' @seealso breedFst USCaucs
#' @importFrom utils write.csv
#' @examples
#' 
#' data(USCaucs)
#' \dontrun{
#'   writeCSV("USCaucs.csv", USCaucs)
#' }
#' 
writeCSV = function(fileName, pop, n = 100, delim = ','){
  numLoci = length(pop$loci)
  Alleles = unique(unlist(lapply(pop$freqs, function(locus){
    as.numeric(names(locus))
  })))
  
  if(any(is.na(Alleles))){
    Alleles[is.na(Alleles)] = 108.1 ## Rares 
  }
  
  Alleles = sort(Alleles)
  numAlleles = length(Alleles)
  
  tbl = matrix(0, nrow = numAlleles, ncol = numLoci)
  colnames(tbl) = pop$loci

  l = 1
  for(loc in pop$freqs){
    a = match(as.numeric(names(loc)), Alleles)
    tbl[a, l] = loc
    l = l + 1
  }
  
  tbl = as.data.frame(rbind(tbl, rep(n, numLoci)))
  tbl$Allele = c(as.character(Alleles), "N")
  tbl = tbl[,c(numLoci + 1, 1:numLoci)]
  write.csv(tbl, file = fileName, row.names = FALSE)
  
  invisible(tbl)
}
