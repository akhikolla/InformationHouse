#' Read PLINK text file.
#'
#' \code{read.plink.text} reads PLINK text files in either the original or transposed format.
#'
#' The PLINK pedigree file should be supplied with the appropriate extension. The corresponding map file can be omitted if it has the same file name as the pedigree file and has the appropriate extension.
#'
#' @param ped,map PLINK files with appropraite extensions.
#' @param output.allele logical, default is to output genotype as alleles.
#' @param na.strings Character vector, set of characters to be treated as missing values.
#' @return A list of three elements: \code{genotype}, \code{fam} and \code{map}.
#' To be consistent with PLINK .ped file, \code{genotype} by default is a n_subject by (2 x n_marker) matrix of alleles, where 1 represents the reference allele and 2 the alternate allele. Alternatively, genotype can be outputted as 0, 1 or 2 copies of reference allele count by using \code{output.allele = FALSE}. Missing values are -9.
#' \code{fam} is a dataframe that contains the first six columns of a PLINK \code{.ped} file.
#' \code{map} is a dataframe that contains the four columns of a PLINK \code{.map} file, with two additional columns: \code{allele_1} for the reference allele type, \code{allele_2} for the alternate allele type.
#' @export
read.plink.text = function(ped, map = NULL, output.allele = TRUE, na.strings = c("0", "-9")){
  nc = nchar(ped)
  ext = substring(ped, nc-3, nc)

  if(is.null(map)){
    if(ext == ".ped"){
      map = paste(substring(ped, 1, nc-3), "map", sep="")
    }else{
      map = paste(substring(ped, 1, nc-4), "tfam", sep="")
    }
  }

  if(ext == ".ped"){
    df.map = read.table(map, na.strings = na.strings, as.is = TRUE)
    if(ncol(df.map) == 4){
      names(df.map) = c("chromosome", "SNP", "cM", "bp")
    }else{
      warning("Non-standard .map file!\n")
    }
    nsnp = nrow(df.map)

    df.ped = read.table(ped, na.strings = na.strings, as.is = TRUE)
    if(ncol(df.ped) == (6 + 2*nsnp)){
      df.fam = df.ped[, 1:6]
      names(df.fam) = c("pedigree", "member", "father", "mother", "sex", "phenotype")
    }else{
      warning("Non-standard .ped file!\n")
    }
    nind = nrow(df.fam)

    if(output.allele){
      L = process_snptext_cpp(df.ped[,7:ncol(df.ped)], 0, 1)
    }else{
      L = process_snptext_cpp(df.ped[,7:ncol(df.ped)], 0, 1, outallele = 0)
    }
  }else{
    df.fam = read.table(map, na.strings = na.strings, as.is = TRUE)
    if(ncol(df.fam) ==  6){
      names(df.fam) = c("pedigree", "member", "father", "mother", "sex", "phenotype")
    }else{
      warning("Non-standard .tfam file!\n")
    }
    nind = nrow(df.fam)

    df.ped = read.table(ped, na.strings = na.strings, as.is = TRUE)
    if(ncol(df.ped) == (4 + 2*nind)){
      df.map = df.ped[,1:4]
      names(df.map) = c("chromosome", "SNP", "cM", "bp")
    }else{
      warning("Non-standard .tped file!\n")
    }
    nsnp = nrow(df.map)

    if(output.allele){
      L = process_snptext_cpp(df.ped[,5:ncol(df.ped)], 1, 1)
    }else{
      L = process_snptext_cpp(df.ped[,5:ncol(df.ped)], 1, 1, outallele = 0)
    }
  }
  output = list(genotype = L[[1]], fam = df.fam, map = cbind(df.map, L[[2]]))
  return(output)
}
