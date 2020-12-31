#' Read PLINK binary file.
#'
#' \code{read.plink.binary} reads PLINK binary \code{.bed} file and the corresponding \code{.bim} and \code{.fam} file.
#'
#' When the three files have the same name, only the \code{.bed} file needs to be specified.
#'
#' @param bed,bim,fam PLINK files with appropriate extensions.
#' @param na.strings string vector, text entries to be treated as NA's.
#' @return A list of three elements: \code{genotype}, \code{fam} and \code{map}.
#' To be consistent with PLINK .bed file, \code{genotype} is a n_subject by n_marker matrix of counts of reference alleles. Missing values are -9.
#' \code{fam} is a dataframe that contains the first six columns of a PLINK \code{.ped} file.
#' \code{map} is a dataframe that contains the four columns of a PLINK \code{.map} file, with two additional columns: \code{allele_1} for the reference allele type, \code{allele_2} for the alternate allele type.
#' @export
read.plink.binary = function(bed, bim = NULL, fam = NULL, na.strings = c("0", "-9")){
  if(is.null(bim)){
    nc = nchar(bed)
    file = substr(bed, 1, nc-4)
    bim = paste(file, ".bim", sep = "")
    if(is.null(fam)){
      fam = paste(file, ".fam", sep = "")
    }
  }

  # read .bim file
  df.bim = read.table(bim, na.strings = na.strings, as.is = TRUE)
  if(ncol(df.bim) == 6){
    names(df.bim) = c("chromosome", "SNP", "cM", "bp", "allele_1", "allele_2")
  }else{
    warning("Non-standard .bim file!")
  }
  nsnp = nrow(df.bim)
  if(length(table(df.bim$SNP)) < nsnp){
    stop("Duplicated SNP name(s)!")
  }

  # read .fam file
  df.fam = read.table(fam, na.strings = na.strings, as.is = TRUE)
  if(ncol(df.fam) == 6){
    names(df.fam) = c("pedigree", "member", "father", "mother", "sex", "phenotype")
  }else{
    warning("Non-standard .fam file!")
  }
  nind = nrow(df.fam)

  geno = read_bed_cpp(bed, nind, nsnp)

  output = list(genotype = geno, fam = df.fam, map = df.bim)
  return(output)
}
