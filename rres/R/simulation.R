#' Simulate inheritance on a given pedigree.
#'
#' \code{sim.recomb} returns inheritance information simulated on a given pedigree over the specified segment length.
#'
#' \code{pedinfo} must contain at least the following components: unique individual ID named \code{member}, father and mother ID named \code{father} and \code{mother}, and sex (1 for male, 2 for female) named \code{sex}. Parents must precede offsprings. Pedigree founders are treated as unrelated.
#'
#' \code{seglength} represents length of genomic segment in Haldane centiMorgan. Recombination breakpoints are simulated under a homogeneous Poisson process with rate \code{seglength}/100.
#'
#' @importFrom stats cor rpois runif
#' @importFrom utils read.table tail
#' @param pedinfo dataframe.
#' @param seglength positive real number.
#' @return A list of matrices for each meiosis. Each matrix has two columns: founder genome labels (fgl) and recombination breakpoints (recomb). Paternal meiosis precedes maternal meiosis.
#' @examples
#' # a simple pedigree with sibling marriage
#' pedigree = as.character(rep(1, 5))
#' member = as.character(c(11, 12, 21, 22, 31))
#' sex = as.numeric(c(1, 2, 1, 2, 1))
#' father = as.character(c(NA, NA, 11, 11, 21))
#' mother = as.character(c(NA, NA, 12, 12, 22))
#' pedinfo = data.frame(pedigree, member, sex, father, mother, stringsAsFactors = FALSE)
#'
#' # simulate inheritance over a segment of 100 centiMorgan
#' sim.recomb(pedinfo, 100)
#' @export
sim.recomb = function(pedinfo, seglength){
  member = father = mother = sex = NULL
  pedinfo = transform(pedinfo, member = as.character(member), father = as.character(father), mother = as.character(mother))

  # check pedigree information
  check.pedinfo(pedinfo)

  meis.count = 1
  output = list()

  for(i in 1:length(pedinfo$member)){
    if(is.na(pedinfo$father[i])){ # founder
      # paternal
      output[[2 * i - 1]] = matrix(c(meis.count, seglength), 1, 2, byrow = TRUE, dimnames = list(NULL, c("fgl", "recomb")))
      # maternal
      output[[2 * i]] = matrix(c(meis.count+1, seglength), 1, 2, byrow = TRUE, dimnames = list(NULL, c("fgl", "recomb")))
      meis.count = meis.count + 2
    }else{ # non-founder
      # paternal inheritance
      fa.id = which(pedinfo$member == pedinfo$father[i]) # father's id
      fa.recomb = sort(runif(rpois(1, seglength/100), min = 0, max = seglength)) # paternal recombination breakpoints
      fa.parental = rep(sample(c(0,1), replace = FALSE), ceiling((length(fa.recomb) + 1)/2))

      if(length(fa.recomb) == 0){
        output[[2 * i - 1]] = output[[fa.id * 2 - fa.parental[1]]]
      }else{
        fa.inheritance = NULL
        fa.switch.count = 1
        fa.gp.count = 1
        fa.gm.count = 1
        fa.inher.count = 1

        while(fa.switch.count < length(fa.recomb) + 2){
          if(fa.parental[fa.switch.count] == 1){ # paternal
            while(fa.switch.count > 1 && output[[fa.id * 2 - 1]][fa.gp.count, 2] < fa.recomb[fa.switch.count - 1]){
              fa.gp.count = fa.gp.count + 1
            }

            if(fa.switch.count > length(fa.recomb)){
              while(fa.gp.count <= nrow(output[[fa.id * 2 - 1]])){
                if(fa.inher.count > 2 && fa.inheritance[fa.inher.count - 2] == output[[fa.id * 2 - 1]][[fa.gp.count]]){
                  fa.inher.count = fa.inher.count - 2
                }
                fa.inheritance[fa.inher.count] = output[[fa.id * 2 - 1]][[fa.gp.count, 1]]
                fa.inheritance[fa.inher.count + 1] = output[[fa.id * 2 - 1]][[fa.gp.count, 2]]
                fa.inher.count = fa.inher.count + 2
                fa.gp.count = fa.gp.count + 1
              }
              break
            }

            while(fa.gp.count <= nrow(output[[fa.id * 2 - 1]])){
              if(fa.inher.count > 2 && fa.inheritance[fa.inher.count - 2] == output[[fa.id * 2 - 1]][[fa.gp.count, 1]]){
                fa.inher.count = fa.inher.count - 2
              }

              if(output[[fa.id * 2 - 1]][fa.gp.count, 2] < fa.recomb[fa.switch.count]){
                fa.inheritance[fa.inher.count] = output[[fa.id * 2 - 1]][[fa.gp.count, 1]]
                fa.inheritance[fa.inher.count + 1] = output[[fa.id * 2 - 1]][[fa.gp.count, 2]]
                fa.inher.count = fa.inher.count + 2
                fa.gp.count = fa.gp.count + 1
              }else{
                fa.inheritance[fa.inher.count] = output[[fa.id * 2 - 1]][[fa.gp.count, 1]]
                fa.inheritance[fa.inher.count + 1] = fa.recomb[fa.switch.count]
                fa.inher.count = fa.inher.count + 2
                break
              }
            }
          }else{ # maternal
            while(fa.switch.count > 1 && output[[fa.id * 2]][fa.gm.count, 2] < fa.recomb[fa.switch.count - 1]){
              fa.gm.count = fa.gm.count + 1
            }

            if(fa.switch.count > length(fa.recomb)){
              while(fa.gm.count <= nrow(output[[fa.id * 2]])){
                if(fa.inher.count > 2 && fa.inheritance[fa.inher.count - 2] == output[[fa.id * 2]][[fa.gm.count]]){
                  fa.inher.count = fa.inher.count - 2
                }
                fa.inheritance[fa.inher.count] = output[[fa.id * 2]][[fa.gm.count, 1]]
                fa.inheritance[fa.inher.count + 1] = output[[fa.id * 2]][[fa.gm.count, 2]]
                fa.inher.count = fa.inher.count + 2
                fa.gm.count = fa.gm.count + 1
              }
              break
            }

            while(fa.gm.count <= nrow(output[[fa.id * 2]])){
              if(fa.inher.count > 2 && fa.inheritance[fa.inher.count - 2] == output[[fa.id * 2]][[fa.gm.count, 1]]){
                fa.inher.count = fa.inher.count - 2
              }

              if(output[[fa.id * 2]][fa.gm.count, 2] < fa.recomb[fa.switch.count]){
                fa.inheritance[fa.inher.count] = output[[fa.id * 2]][[fa.gm.count, 1]]
                fa.inheritance[fa.inher.count + 1] = output[[fa.id * 2]][[fa.gm.count, 2]]
                fa.inher.count = fa.inher.count + 2
                fa.gm.count = fa.gm.count + 1
              }else{
                fa.inheritance[fa.inher.count] = output[[fa.id * 2]][[fa.gm.count, 1]]
                fa.inheritance[fa.inher.count + 1] = fa.recomb[fa.switch.count]
                fa.inher.count = fa.inher.count + 2
                break
              }
            }
          }
          fa.switch.count = fa.switch.count + 1
        }
        output[[2 * i - 1]] = matrix(fa.inheritance, ncol = 2, byrow = TRUE, dimnames = list(NULL, c("fgl", "recomb")))
      }

      # maternal inheritance
      mo.id = which(pedinfo$member == pedinfo$mother[i]) # mother's id
      mo.recomb = sort(runif(rpois(1, seglength/100), min = 0, max = seglength)) # paternal recombination breakpoints
      mo.parental = rep(sample(c(0,1), replace = FALSE), ceiling((length(mo.recomb) + 1)/2))

      if(length(mo.recomb) == 0){
        output[[2 * i]] = output[[mo.id * 2 - mo.parental[1]]]
      }else{
        mo.inheritance = NULL
        mo.switch.count = 1
        mo.gp.count = 1
        mo.gm.count = 1
        mo.inher.count = 1

        while(mo.switch.count < length(mo.recomb) + 2){
          if(mo.parental[mo.switch.count] == 1){ # paternal
            while(mo.switch.count > 1 && output[[mo.id * 2 - 1]][mo.gp.count, 2] < mo.recomb[mo.switch.count - 1]){
              mo.gp.count = mo.gp.count + 1
            }

            if(mo.switch.count > length(mo.recomb)){
              while(mo.gp.count <= nrow(output[[mo.id * 2 - 1]])){
                if(mo.inher.count > 2 && mo.inheritance[mo.inher.count - 2] == output[[mo.id * 2 - 1]][[mo.gp.count]]){
                  mo.inher.count = mo.inher.count - 2
                }
                mo.inheritance[mo.inher.count] = output[[mo.id * 2 - 1]][[mo.gp.count, 1]]
                mo.inheritance[mo.inher.count + 1] = output[[mo.id * 2 - 1]][[mo.gp.count, 2]]
                mo.inher.count = mo.inher.count + 2
                mo.gp.count = mo.gp.count + 1
              }
              break
            }

            while(mo.gp.count <= nrow(output[[mo.id * 2 - 1]])){
              if(mo.inher.count > 2 && mo.inheritance[mo.inher.count - 2] == output[[mo.id * 2 - 1]][[mo.gp.count, 1]]){
                mo.inher.count = mo.inher.count - 2
              }

              if(output[[mo.id * 2 - 1]][mo.gp.count, 2] < mo.recomb[mo.switch.count]){
                mo.inheritance[mo.inher.count] = output[[mo.id * 2 - 1]][[mo.gp.count, 1]]
                mo.inheritance[mo.inher.count + 1] = output[[mo.id * 2 - 1]][[mo.gp.count, 2]]
                mo.inher.count = mo.inher.count + 2
                mo.gp.count = mo.gp.count + 1
              }else{
                mo.inheritance[mo.inher.count] = output[[mo.id * 2 - 1]][[mo.gp.count, 1]]
                mo.inheritance[mo.inher.count + 1] = mo.recomb[mo.switch.count]
                mo.inher.count = mo.inher.count + 2
                break
              }
            }
          }else{ # maternal
            while(mo.switch.count > 1 && output[[mo.id * 2]][mo.gm.count, 2] < mo.recomb[mo.switch.count - 1]){
              mo.gm.count = mo.gm.count + 1
            }

            if(mo.switch.count > length(mo.recomb)){
              while(mo.gm.count <= nrow(output[[mo.id * 2]])){
                if(mo.inher.count > 2 && mo.inheritance[mo.inher.count - 2] == output[[mo.id * 2]][[mo.gm.count]]){
                  mo.inher.count = mo.inher.count - 2
                }
                mo.inheritance[mo.inher.count] = output[[mo.id * 2]][[mo.gm.count, 1]]
                mo.inheritance[mo.inher.count + 1] = output[[mo.id * 2]][[mo.gm.count, 2]]
                mo.inher.count = mo.inher.count + 2
                mo.gm.count = mo.gm.count + 1
              }
              break
            }

            while(mo.gm.count <= nrow(output[[mo.id * 2]])){
              if(mo.inher.count > 2 && mo.inheritance[mo.inher.count - 2] == output[[mo.id * 2]][[mo.gm.count, 1]]){
                mo.inher.count = mo.inher.count - 2
              }

              if(output[[mo.id * 2]][mo.gm.count, 2] < mo.recomb[mo.switch.count]){
                mo.inheritance[mo.inher.count] = output[[mo.id * 2]][[mo.gm.count, 1]]
                mo.inheritance[mo.inher.count + 1] = output[[mo.id * 2]][[mo.gm.count, 2]]
                mo.inher.count = mo.inher.count + 2
                mo.gm.count = mo.gm.count + 1
              }else{
                mo.inheritance[mo.inher.count] = output[[mo.id * 2]][[mo.gm.count, 1]]
                mo.inheritance[mo.inher.count + 1] = mo.recomb[mo.switch.count]
                mo.inher.count = mo.inher.count + 2
                break
              }
            }
          }
          mo.switch.count = mo.switch.count + 1
        }
        output[[2 * i]] = matrix(mo.inheritance, ncol = 2, byrow = TRUE, dimnames = list(NULL, c("fgl", "recomb")))
      }
    }
  }
  return(output)
}


#' Simulate artificial haplotypes.
#'
#' \code{sim.haplotype} returns haplotypes of the specified number of SNPs simulated under linkage equilibrium.
#'
#' \code{freq} are reference allele frequencies. \code{nhaplo} haplotypes are simulated independently.
#'
#' @param freq vector of values between 0 and 1.
#' @param nhaplo positive integer.
#' @return A matrix of \code{nhaplo} rows and \code{length(freq)} columns. Reference alleles are coded 1, alternate alleles are coded 2.
#' @examples
#' nsnp = 7 # number of SNPs
#' freq = runif(nsnp, 0.05, 0.95)
#' nhaplo = 4 # number of founder haplotypes
#' sim.haplotype(freq, nhaplo)
#' @export
sim.haplotype = function(freq, nhaplo){
  output = sapply(freq, function(x) sample(c(1, 2), nhaplo, prob = c(x, 1-x), replace = TRUE))
  return(output)
}


#' Populate SNPs.
#'
#' \code{populate.snp} assigns alleles to markers, given inheritance information and founder haplotypes.
#'
#' \code{inheritance} is a list of matrices produced by, e.g., \code{\link{sim.recomb}}. Each matrix contains a column of founder genome labels and a column of recombination breakpoints for the corresponding meiosis.
#'
#' \code{haplotype} is a numeric matrix. The matrix is number of haplotypes by number of markers in dimension. Standard coding in this package is 1 for reference allele and 2 for alternate allele. This coding is required when \code{output.allele = FALSE}. Input data with different coding of alleles can be recoded using \code{\link{recode.snpdata}}. Number of haplotypes cannot be fewer than the number of founder genome labels in \code{inheritance}. The haplotypes will be assigned to each founder genome label in given order.
#'
#' \code{marker} is a vector of marker genetic positions in Haldane centiMorgan in ascending order. Range of marker positions cannot exceed range covered by \code{inheritance}.
#'
#' \code{member.index} contains indices of members in the pedigree that we wish to output data. Default value is \code{FALSE}, in which case marker data on everyone will be produced. \code{\link{get.pedindex}} can help find indices given member ID.
#'
#' \code{output.allele} determines if one or two numbers will be used to represent data at each marker. Default is \code{TRUE}, in which case marker data is represented by two ordered (paternal first) alleles. Otherwise marker data is represented by a single number (0, 1 or 2) of reference alleles.
#'
#' \code{output.haplotype} determines if haplotype data are separate in output. It is only used when \code{output.allele = TRUE}. Default value is \code{FALSE}, in which case each row in the output matrix represents ordered genotypes from all markers of the same individual. Otherwise each row in the output matrix represents a parental haplotype.
#'
#' @param inheritance list of numeric matrices.
#' @param haplotype numeric matrix.
#' @param marker numeric vector.
#' @param member.index integer vector.
#' @param output.allele logical.
#' @param output.haplotype logical.
#' @return A matrix of genotypic/haplotypic data. The matrix is in individual major, where marker data for each individual/meiosis are found on the same row. Exact format of the matrix depends on various input arguments.
#' @examples
#' # a simple pedigree with sibling marriage
#' pedigree = as.character(rep(1, 5))
#' member = as.character(c(11, 12, 21, 22, 31))
#' sex = as.numeric(c(1, 2, 1, 2, 1))
#' father = as.character(c(NA, NA, 11, 11, 21))
#' mother = as.character(c(NA, NA, 12, 12, 22))
#' pedinfo = data.frame(pedigree, member, sex, father, mother, stringsAsFactors = FALSE)
#'
#' L = 100.0 # segment length
#' nsnp = 10 # number of SNPs
#' nhaplo = 4 # number of founder haplotypes
#' inher = sim.recomb(pedinfo, L)
#' haplo = matrix(c(3,4,4,4), nhaplo, nsnp)
#' marker = sort(runif(nsnp, 0, L))
#'
#' # output genotype data for the 4th and 5th member
#' # of pedigree, genotype data displayed as two alleles
#' populate.snp(inher, haplo, marker, c(4, 5))
#'
#' # output haplotype data for the 4th and 5th member of pedigree
#' populate.snp(inher, haplo, marker,c(4, 5), output.haplotype = TRUE)
#'
#' # output genotype data for all members, genotype data
#' # displayed as counts of reference alleles
#' geno = recode.snpdata(haplo, input.haplotype = TRUE, output.haplotype = TRUE)[[1]]
#' populate.snp(inher, geno, marker, output.allele = FALSE)
#' @export
#' @useDynLib rres
#' @importFrom Rcpp sourceCpp
populate.snp = function(inheritance, haplotype, marker, member.index = NULL, output.allele = TRUE, output.haplotype = FALSE){
  nhaplo = length(inheritance) # total number of haplotypes
  nsnp = length(marker) # total number of snps on genomic segment

  # get number of founder haplotypes
  fgl = NULL
  nfounder2 = length(unique(lapply(inheritance, function(x) c(fgl, x[1,1])))) # number of founder haplotypes

  # get indices of haplotypes that we wish to output
  if(!is.null(member.index)){
    haplo.index = NULL
    for(i in member.index){
      haplo.index = c(haplo.index, c(i*2-1, i*2))
    }
  }else{
    haplo.index = c(1:length(inheritance))
  }
  nind = length(haplo.index)/2

  # check if enough founder haplotypes were supplied
  if(nrow(haplotype) < nfounder2){
    stop("Insufficient number of founder haplotypes.")
  }

  L = inheritance[[1]][1,2]
  if(min(marker) < 0 || max(marker) > L){
    stop("Marker positon out of bounds.")
  }

  if(ncol(haplotype) != nsnp){
    stop("Haplotype input and marker position input contain different number of markers.")
  }

  haplo.mat = matrix(0, nind*2, nsnp)
  for(i in 1:length(haplo.index)){
    current.index = 1
    current.fgl = inheritance[[haplo.index[i]]][current.index, 1]
    current.recomb = inheritance[[haplo.index[i]]][current.index, 2]
    start.mindex = 1
    end.mindex = 1

    if(current.recomb == L){
      haplo.mat[i, ] = haplotype[current.fgl, ]
      next
    }

    while(current.recomb < L && end.mindex <= nsnp){
      if(marker[end.mindex] > current.recomb){
        if(start.mindex != end.mindex){
          haplo.mat[i, start.mindex:(end.mindex-1)] = haplotype[current.fgl, start.mindex:(end.mindex-1)]
          start.mindex = end.mindex
        }
        current.index = current.index + 1
        current.fgl = inheritance[[haplo.index[i]]][current.index, 1]
        current.recomb = inheritance[[haplo.index[i]]][current.index, 2]
      }else{
        end.mindex = end.mindex + 1
      }
    }
    haplo.mat[i, start.mindex:nsnp] = haplotype[current.fgl, start.mindex:nsnp]
  }

  if(output.allele){
    if(output.haplotype){
      return(haplo.mat)
    }else{
      geno.mat = t(apply(matrix(as.vector(t(haplo.mat)), nind, nsnp*2, byrow = TRUE), 1, function(x) as.vector(matrix(x, nrow = 2, byrow = TRUE))))
      return(geno.mat)
    }
  }else{
    geno.mat = t(sapply(c(1:nind), function(x) 4 - haplo.mat[2*x-1,] - haplo.mat[2*x,]))
    return(geno.mat)
  }
}
