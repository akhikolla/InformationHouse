#' CODIS STR Loci allele frequency data
#' 
#' This data structure
#' 
#' 
#' @name fbiCaucs
#' @docType data
#' @format This data set is a list which has two sub-lists. The lists are named
#' loci and freqs. loci is a vector of the 13 CODIS STR locus names. freqs is a
#' list of 13 vectors, each vector contains the allele frequencies published
#' for US Caucasians in Budowle et al. (2001).
#' @author James M. Curran
#' @seealso USCaucs
#' @references Budowle B, Shea B, Niezgoda S, Chakraborty R. (2001),
#' \emph{CODIS STR loci data from 41 sample populations}, J. Forensic Sci.
#' 46:453-89.
#' @keywords datasets
#' @examples
#' 
#' data(fbiCaucs)
#' names(fbiCaucs)
#' fbiCaucs$loci
#' names(fbiCaucs$freqs)
#' fbiCaucs$freqs[[1]]
#' names(fbiCaucs$freqs[[1]])
#' fbiCaucs$freqs[[1]][1]
#' 
NULL

#' CODIS STR Loci allele frequency data
#' 
#' This data structure
#' 
#' 
#' @name USCaucs
#' @docType data
#' @format This data set is a list which has two sub-lists. The lists are named
#' loci and freqs. loci is a vector of the 13 CODIS STR locus names. freqs is a
#' list of 13 vectors, each vector contains the allele frequencies published
#' for US Caucasians in Budowle and Moretti (1999).  The raw data is available
#' from
#' \url{http://www.fbi.gov/about-us/lab/forensic-science-communications/fsc/july1999/dnaloci.txt}
#' @author James M. Curran
#' @seealso fbiCaucs
#' @references Budowle, B. and Moretti, T.R. (1999), \emph{Genotype Profiles
#' for Six Population Groups at the 13 CODIS Short Tandem Repeat Core Loci and
#' Other PCR Based Loci}, Forensic Science Communications 1(2).
#' @keywords datasets
#' @examples
#' 
#' data(USCaucs)
#' names(USCaucs)
#' USCaucs$loci
#' names(USCaucs$freqs)
#' USCaucs$freqs[[1]]
#' names(USCaucs$freqs[[1]])
#' USCaucs$freqs[[1]][1]
NULL