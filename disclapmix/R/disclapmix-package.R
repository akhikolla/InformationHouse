.onUnload <- function(libpath) {
  library.dynam.unload("disclapmix", libpath)
}

#' Y-STR haplotypes
#' 
#' 185 Y-STR 10 loci haplotypes
#' 
#' 
#' @name danes
#' @docType data
#' @format A data frame with 185 observations on the following 10 loci
#' (\code{n} is the number of times each haplotype has been observed)
#' \describe{ \item{DYS19}{} \item{DYS389I}{}
#' \item{DYS389II}{} \item{DYS390}{} \item{DYS391}{}
#' \item{DYS392}{} \item{DYS393}{} \item{DYS437}{}
#' \item{DYS438}{} \item{DYS439}{} \item{n}{} }
#' @source ''Y-chromosome STR haplotypes Danes'' by Hallenberg et al (2005),
#' http://www.sciencedirect.com/science/article/pii/S0379073805000095
#' @keywords datasets
NULL


#' disclapmix
#' 
#' Discrete Laplace Mixture Inference using the EM Algorithm
#' 
#' @docType package
#' @author Mikkel Meyer Andersen <mikl@math.aau.dk> 
#' @import Rcpp disclap
#' @importFrom Rcpp evalCpp
#' @importFrom cluster pam clara
#' @importFrom MASS mvrnorm
#' @importFrom stats simulate predict as.dist glm.control glm.fit hclust median sd model.matrix
#' @importFrom graphics plot
#' @importFrom methods is
#' @importFrom utils head tail
#' @useDynLib disclapmix
#' @name disclapmix
NULL

