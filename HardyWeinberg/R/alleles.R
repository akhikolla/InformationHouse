alleles <- function (x,fromlabels=TRUE) 
{
  if(fromlabels) {
    geno.unique <- sort(unique(names(x)))  
  } else {
    geno.unique <- sort(unique(x))
  }
  hemizygous.genotypes <- geno.unique[nchar(geno.unique)==1]
  diploid.genotypes <- geno.unique[nchar(geno.unique)==2]
  diploid.alleles <- sort(unique(c(substr(diploid.genotypes, 1, 1), 
                                     substr(diploid.genotypes, 2, 2))))
  alleles <- sort(union(hemizygous.genotypes,diploid.alleles))
  if(any(nchar(geno.unique)>2)) stop("alleles: there are genotypes with more that 2 constituent alleles")
  if(any(nchar(geno.unique)==0)) stop("alleles: there are genotypes with empty alleles")
  return(alleles)
}
