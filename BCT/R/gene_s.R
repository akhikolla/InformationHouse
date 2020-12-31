#' SARS-CoV-2 gene S
#' 
#' This dataset contains the spike (S) gene, in positions 21,563–25,384 of the SARS-CoV-2 genome. 
#' The importance of this gene is that it codes for the surface glycoprotein whose function was identified in Yan et al. (2020) and Lan et al. (2020) as critical, 
#' in that it binds onto the Angiotensin Converting Enzyme 2 (ACE2) receptor on human epithelial cells, 
#' giving the virus access to the cell and thus facilitating the COVID-19 disease. The gene sequence is mapped to the alphabet {0,1,2,3} via the obvious map A->0, C->1, G->2, T->3.
#' @docType data
#' @format An object of class \code{"character"}. 
#'
#' @keywords datasets
#'
#' @references {Wu, F., S. Zhao, B. Yu, et al. (2020). A new coronavirus associated with human respiratory disease in China. Nature 579(7798), 265–269.}
#' (\href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7094943/}{PMC})
#' @references {Yan, R., Y. Zhang, Y. Li, L. Xia, Y. Guo, and Q. Zhou (2020). Structural basis for the recognition of SARS-CoV-2 by full-length human ACE2. Science 367(6485), 1444–1448.}
#' (\href{https://science.sciencemag.org/content/367/6485/1444}{aaas})
#' @references {Lan, J., J. Ge, J. Yu, et al. (2020). Structure of the SARS-CoV-2 spike receptor-binding domain bound to the ACE2 receptor. Nature 581, 215–220.}
#' (\href{https://www.nature.com/articles/s41586-020-2180-5}{nature})
#' @examples
#' BCT(gene_s, 5)

"gene_s"