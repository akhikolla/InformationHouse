#' Read a set of profiles from a file
#'
#' Reads a set of profiles from a file
#'
#' The alleles are recorded integers rather than the STR repeat numbers. This
#' speeds up computation immensely when calculating any of the LRs or IBS.
#'
#' @param fileName a path to the profile file.
#' @param freqs A list containing two lists labelled loci and freqs. The second
#' list is a list of vectors containing the allele frequencies of each allele
#' at each locus in the multiplex. If this is left \code{NULL}, then it is calculated from 
#' the profile file.
#' @param sep a character that delimits the fields in the profile file.
#' @param header a boolean witch is \code{TRUE} if the profile file has a column header line.
#' @param id a column number indicating which column the profile id's are stored. If \code{id == -1}, then this means there is no id information.
#' @param discardMissing if \code{TRUE}, then all profiles which have alleles which cannot be matched to the frequency file are ommitted and returned
#' in the return list.
#' @return a list containing a \code{data.frame} of profiles where the alleles have been recoded 
#' to the allele index number, rather than the allele itself, and a set of frequencies in the same
#' format as you would get from \code{\link{readFreqs}}. If \code{freqs} have been supplied, then this 
#' will just be the same set of frequencies, if they have not, then this will be calculated
#' from the profiles. Given that the profiles generally do not have any locus name information
#' the loci will just be labelled Locus1, Locus2, \ldots. If there are missing values then the raw missing profiles are returned
#' @author James M. Curran
#'
#' @importFrom utils read.table
#' @importFrom stats complete.cases
#' @export readProfiles
readProfiles = function(fileName,
                        freqs = NULL,
                        sep = "\t",
                        header = FALSE,
                        id = 1,
                        discardMissing = TRUE) {
  if (header) {
    rawProfiles = read.table(fileName,
                      sep = sep,
                      header = FALSE,
                      skip = 1)
  } else{
    rawProfiles = read.table(fileName, sep = sep, header = FALSE)
  }
  
  cat(paste0("Read ", nrow(rawProfiles), " profiles from ", fileName, "\n"))
  
  prof = rawProfiles
  
  idCol = NULL
  if (id != -1) {
    idCol = rawProfiles[, id]
    prof = prof[, -id]
  }
  
  if(is.null(freqs)){
    cat("No frequency file supplied so compiling from raw profiles\n")
    
    if(ncol(rawProfiles) %% 2 != 0){
      stop("The profile file must contain an even number of columns excluding the id column\n")
    }else{
      numLoci = ncol(rawProfiles) / 2
      
      freqs$loci = paste("Locus", 1:numLoci, sep = "_")
      freqs$freqs = vector(mode = "list", length = numLoci)
      names(freqs$freqs) = freqs$loci
      
      for(nloc in 1:numLoci){
        i1 = 2 * nloc - 1
        i2 = i1 + 1
        freqs$freqs[[nloc]] = table(c(rawProfiles[,i1], rawProfiles[,i2]))
        freqs$freqs[[nloc]] = freqs$freqs[[nloc]] / sum(freqs$freqs[[nloc]])
      }
    }
  }
  
  numLoci = length(freqs$freqs)
  
  for (loc in 1:numLoci) {
    locus = freqs[["freqs"]][[loc]]
    Alleles = as.numeric(names(locus))
    
    i1 = 2 * loc - 1
    i2 = i1 + 1
    prof[, i1] = match(prof[, i1], Alleles)
    prof[, i2] = match(prof[, i2], Alleles)
  }
  
  rownames(prof) = 1:nrow(prof)
  missing = !complete.cases(prof)
  if(any(missing)){
    cat(paste0("File contained ", sum(missing), " profiles with alleles that couldn't be resovled from frequency file\n"))
    cat("These profiles will be ommitted from analysis, but are returned in is object as $missing\n")
    return(list(profiles = prof[!missing, ], freqs = freqs, missing = list(raw = rawProfiles[missing,], mapped = prof[missing, ])))
  }
  
  return(list(profiles = rawProfiles, freqs = freqs))
}
