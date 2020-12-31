#' SARS-CoV-2 genome
#' 
#' The severe acute respiratory syndrome coronavirus, SARS-CoV-2, is the novel coronavirus responsible for the Covid-19 global pandemic in 2019-20. This dataset contains the SARS-CoV-2 genome, 
#' available in the GenBank database as the sequence MN908947.3.
#' It consists of n = 29903 base pairs. The gene sequence is mapped to the alphabet {0,1,2,3} via the obvious map A->0, C->1, G->2, T->3.
#'  
#' @docType data
#' @format An object of class \code{"character"}. 
#'
#' @keywords datasets
#'
#' @references {K. Clark, I. Karsch-Mizrachi, D.J. Lipman, J. Ostell, and E.W. Sayers. GenBank. Nucleic Acids Research, 44(D1): D67–D72, January 2016.}
#' (\href{https://www.ncbi.nlm.nih.gov}{ncbi})
#' @references {F. Wu, S. Zhao, B. Yu, et al. A new coronavirus associated with human respiratory disease in China. Nature, 579(7798):265–269, Februrary 2020.}
#' (\href{https://pubmed.ncbi.nlm.nih.gov/32296181/}{PubMed})
#' @examples
#' BCT(sars_cov_2, 5)

"sars_cov_2"