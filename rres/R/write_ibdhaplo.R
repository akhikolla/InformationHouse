#' Write IBDHAPLO
#' 
#' \code{write.ibdhaplo} prepares the marker data file for running IBDHAPLO.
#' 
#' The input marker data needs to be subject/haplotype by marker/allele. For example, suppose \code{data} is a 4x10 matrix, use \code{input.allele = FALSE} if \code{data} contains counts of reference alleles of 4 individuals at 10 markers; use \code{input.haplotype = TRUE} if \code{data} contains allelic types of 4 haplotypes at 10 markers; use default options if \code{data} contains allelic types of 4 individuals at 5 markers.
#' 
#' @param marker numeric vector, marker genetic positions in cM.
#' @param freq numeric vector, marker reference allele frequencies.
#' @param data numeric matrix, genetic marker data.
#' @param member string vector, member ID.
#' @param input.allele logical, default TRUE.
#' @param input.haplotype logical, default FALSE.
#' @param outfile string, output file name.
#' @references MORGAN Tutorial, \url{https://www.stat.washington.edu/thompson/Genepi/MORGAN/Morgan.shtml}.
#' @references Brown et al. (2012) Genetics 190:1447-1460, \url{https://www.ncbi.nlm.nih.gov/pubmed/22298700}.
#' 
#' @examples 
#' \dontrun{
#' nsnp = 7 # number of SNPs
#' freq = runif(nsnp, 0.05, 0.95)
#' nhaplo = 4 # number of founder haplotypes
#' haplotype = sim.haplotype(freq, nhaplo)
#' marker = sort(runif(7,0,100))
#' write.ibdhaplo(marker, freq, haplotype, member = c("ind1", "ind2"), 
#' input.haplotype = TRUE)
#' }
#' @export
write.ibdhaplo = function(marker, freq, data, member, input.allele = TRUE, input.haplotype = FALSE, outfile = tempfile("ibdhaplo", fileext = ".txt")){

  if(!input.allele){
    genotype = t(apply(data, 1, function(x) c(sapply(x, function(y) c(rep(1,y),rep(2,2-y))))))
  }else{
    if(input.haplotype){
      genotype = recode.snpdata(data, input.haplotype = TRUE, output.haplotype = FALSE)[[1]]
    }
  }
  
  if(length(marker) != length(freq) || length(marker) != ncol(genotype)/2){
    stop("Inconsistent number of SNPs.")
  }
  
  if(nrow(genotype) != length(member)){
    stop("Inconsistent number of individuals.")
  }
  
  write_markers_cpp(marker, freq, genotype, member, outfile)
}



