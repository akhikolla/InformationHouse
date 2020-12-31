#' GRM for a pair of individuals.
#' 
#' \code{grm.pair} computes relatedness estimates between two individuals.
#' 
#' \code{geno1} and \code{geno2} are vectors of counts of reference alleles. \code{freq} is the vector of reference allele frequencies.
#' 
#' The default \code{method} is "twostep", other options include "classic", "robust" and "general". When using the default "twostep" method, user can supply an initial estimate through \code{init.est} to bypass the first step. When "general" is selected, \code{weights} must also be specified. The difference between the two-step GRM, classic GRM and robust GRM is discussed in Wang et al.  (2017).
#' 
#' @references Wang et al. (2017) Genetics 205:1063-1078, \url{https://www.ncbi.nlm.nih.gov/pubmed/28100587}.
#' 
#' @param geno1,geno2 numeric vector.
#' @param freq numeric vector, values between 0 and 1.
#' @param method string.
#' @param weights numeric vector, values between 0 and 1.
#' @param init.est numeric.
#' @return An estimate of realized relatedness.
#' @examples 
#' # simulate genotypes for a full sib pair
#' pedigree = as.character(rep(1, 4))
#' member = as.character(c(11, 12, 21, 22))
#' sex = as.numeric(c(1, 2, 1, 2))
#' father = as.character(c(NA, NA, 11, 11))
#' mother = as.character(c(NA, NA, 12, 12))
#' pedinfo = data.frame(pedigree, member, sex, father, mother, stringsAsFactors = FALSE)
#' set.seed(1)
#' inher = sim.recomb(pedinfo, 3500) # on a hypothetical chromosome
#' 
#' nsnp = 100000
#' marker = seq(0,3500,length.out=nsnp)
#' freq = runif(nsnp, 0.05, 0.95)
#' haplo = sim.haplotype(freq, 4)
#' geno = populate.snp(inher, haplo, marker, output.allele = FALSE)
#' 
#' # simulation truth
#' ibd.proportion(inher,3,4)
#' 
#' # different GRM estimates
#' grm.pair(geno[3,], geno[4,], freq, method = "twostep")
#' grm.pair(geno[3,], geno[4,], freq, method = "classic")
#' grm.pair(geno[3,], geno[4,], freq, method = "robust")
#' grm.pair(geno[3,], geno[4,], freq, method = "general", weights = sample(freq, nsnp)/sum(freq))
#' 
#' # compute the relatedness matrix
#' grm.matrix(geno, freq)
#' grm.matrix(geno, freq, method = "robust")
#' @seealso \code{\link{grm.matrix}}
#' @export
grm.pair = function(geno1, geno2, freq, method = "twostep", weights = NULL, init.est = NULL){
  # check length of input vectors
  if(length(geno1) != length(geno2) || length(geno1) != length(freq)){
    stop("Check input data dimensions.")
  }

  if(method == "classic"){
    est = mean((geno1 - 2 * freq) * (geno2 - 2 * freq) / (2 * freq * (1 - freq)))
  }else if(method == "robust"){
    est = sum((geno1 - 2 * freq) * (geno2 - 2 * freq)) / sum(2 * freq * (1 - freq))
  }else if(method == "twostep"){
    if(is.null(init.est)){
      init.est = mean((geno1 - 2 * freq) * (geno2 - 2 * freq) / (2 * freq * (1 - freq)))
    }
    al = 2 * freq / (1 + init.est) + init.est / (1 + init.est)
    vl = ((al - 2 * freq)^2 + init.est * (al - 1)^2 - init.est/2) / (4 * freq * (1 - freq)) + init.est * (1 - init.est/2)/2 + 1/4
    wl = 1 / vl / (sum(1 / vl))
    est = sum(wl * (geno1 * geno2 - al * (geno1 + geno2) + 4 * al * freq - 4 * freq^2) / (2 * freq * (1 - freq)))
  }else if(method == "general"){
    if(is.null(weights)){
      stop("Weights required.")
    }else{
      if(length(freq) != length(weights)){
        stop("Check length of weights vector.")
      }
      
      if(abs(sum(weights)-1) > 0.000001){
        stop("Weights do not add up to 1.")
      }
    }
    
    est = sum((geno1 - 2 * freq) * (geno2 - 2 * freq) / (2 * freq * (1 - freq)) * weights)
  }else{
    stop("Invalid method.")
  }
  return(est)
}


#' GRM for multiple individuals
#' 
#' \code{grm.matrix} computes relatedness estimates between every pairs of individuals.
#' 
#' \code{genotype} is the matrix of counts of reference alleles. Rows represents subjects and columns represents SNP markers. \code{freq} is the vector of reference allele frequencies.
#' 
#' The default \code{method} is "twostep", other options include "classic", "robust" and "general". When using the default "twostep" method, user can supply an initial estimate through \code{init.est} to bypass the first step. When "general" is selected, \code{weights} must also be specified. The difference between the two-step GRM, classic GRM and robust GRM is discussed in Wang et al.  (2017).
#' 
#' @references Wang et al. (2017) Genetics 205:1063-1078, \url{https://www.ncbi.nlm.nih.gov/pubmed/28100587}.
#' 
#' @param genotype numeric matrix.
#' @param freq numeric vector, values between 0 and 1.
#' @param method string.
#' @param weights numeric vector, values between 0 and 1.
#' @param init.est numeric.
#' @seealso  \code{\link{grm.pair}}.
#' @export
grm.matrix = function(genotype, freq, method = "twostep", weights = NULL, init.est = NULL){
  n = nrow(genotype)
  p = ncol(genotype)
  
  # check input data dimension
  if(n < 2 || p != length(freq)){
    stop("Check input data dimensions.")
  }
  
  if(method == "classic"){
    est = crossprod(apply(genotype, 1, function(x) (x - 2*freq)/sqrt(2*freq*(1-freq))))/p
  }else if(method == "robust"){
    est = tcrossprod((genotype - 2*matrix(freq, n, p, byrow = TRUE)))/2/sum(freq*(1-freq))
  }else if(method == "twostep"){
    if(is.null(init.est)){
      init.est = crossprod(apply(genotype, 1, function(x) (x - 2*freq)/sqrt(2*freq*(1-freq))))/p
    }else{
      if(nrow(init.est) != n || ncol(init.est) != n){
        stop("Check initial estimates matrix dimension.")
      }
    }
    est = twostep_grm_cpp(genotype, freq, init.est)
  }else if(method == "general"){
    if(is.null(weights)){
      stop("Weights required.")
    }else{
      if(length(freq) != length(weights)){
        stop("Check length of weights vector.")
      }
      
      if(abs(sum(weights)-1) > 0.000001){
        stop("Weights do not add up to 1.")
      }
    }
    
    est = tcrossprod((genotype - 2*matrix(freq, n, p, byrow = TRUE))%*%diag(1/sqrt(2*freq*(1-freq)))%*%diag(sqrt(weights)))
  }else{
    stop("Invalid method.")
  }
  return(est)
}

