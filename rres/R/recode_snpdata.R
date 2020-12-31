#' Recode SNP marker data.
#' 
#' \code{recode.snpdata} recodes SNP marker data for use with other functions in this package.
#' 
#' The standard marker data used by other functions of this package takes one of three forms: (a) subjects by row, counts of reference alleles by column; (b) subjects by row, allelic types (2 per marker) by column; (c) haplotypes (2 per subject) by row, allelic types by column. Reference alleles are coded 1, alternate alleles are coded 2.
#' 
#' By default, \code{snp.major = FALSE}, set it to \code{TRUE} if input matrix has SNPs by row and allelic types (2 per subject) by column. \code{ma.ref = FALSE}, set it to \code{TRUE} if the minor allele is to be the reference allele. \code{input.haplotype = FALSE}, set it to \code{TRUE} if input matrix has haplotypes (2 per subject) by row and allelic types by column. \code{output.allele = TRUE}, set it to \code{FALSE} if counts of reference alleles is the desired output format. \code{output.haplotype = FALSE}, set it to \code{TRUE} if recoded marker data by haplotype is the desired output format. \code{input.haplotype} is only invoked when \code{snp.major = FALSE}. \code{output.haplotype} is only invoked when \code{output.allele = TRUE}. 
#' 
#' @param data numeric matrix or dataframe.
#' @param snp.major logical.
#' @param ma.ref logical.
#' @param input.haplotype logical.
#' @param output.allele logical.
#' @param output.haplotype logical.
#' @param na.string numeric or character vector.
#' @return A list of two elements. First element named \code{data} is a matrix of recoded marker data in specified format. Second element is a dataframe named \code{alleles} that specifies reference/alternate alleles at all markers.
#' @examples 
#' test.dat = matrix(c(3,4,4,3), 4, 10)
#' 
#' # treat test.dat as 4 input haplotypes of two subjects at 10 SNP markers,
#' # output recoded data as haplotypes
#' recode.snpdata(test.dat, input.haplotype = TRUE, output.haplotype = TRUE)
#' 
#' # treat test.dat as 4 input haplotypes of two subjects at 10 SNP markers,
#' # output recoded data as counts of reference alleles
#' recode.snpdata(test.dat, input.haplotype = TRUE, output.allele = FALSE)
#' #'
#' # treat test.dat as allelic types at 5 SNPs of 4 subjects,
#' # output recoded data as haplotypes
#' recode.snpdata(test.dat, output.haplotype = TRUE)
#' @export
recode.snpdata = function(data, snp.major = FALSE, ma.ref = FALSE, input.haplotype = FALSE, output.allele = TRUE, output.haplotype = FALSE, na.string = NULL){
  n = nrow(data)
  p = ncol(data)
  
  if(snp.major){
    if(p%%2 == 1){
      stop("Check SNP data dimension.")
    }
  }else{
    if(input.haplotype){
      if(n%%2 == 1){
        stop("Check SNP data dimension.")
      }
    }else{
      if(p%%2 == 1){
        stop("Check SNP data dimension.")
      }
    }
  }
  
  data = as.data.frame(data, stringsAsFactors = FALSE)
  
  if(!is.null(na.string)){
    for(i in na.string){
      data[data == i] = NA
    }
  }
  
  output = process_snptext_cpp(data, sum(snp.major), sum(ma.ref), sum(input.haplotype), sum(output.allele), sum(output.haplotype))
  return(output)
}

