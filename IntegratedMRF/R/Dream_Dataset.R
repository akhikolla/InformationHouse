#' NCI-Dream Drug Sensitivity Prediction Challenge Dataset
#'
#' A demo dataset of different genomic characterizations and drug sensitivity selected for NCI-Dream
#' Drug Sensitivity Prediction Challenge dataset.
#'
#' @format A list of 6 variables containing genomic characterizations and drug sensitivity:
#' \describe{
#'   \item{finalX_Dream}{List of 5 Matrices where the matrices represent different 
#'   genomic characterizations of Gene Expression, Methylation, RNA sequencing, Reverse Phase Protein Array (RPPA) and 
#'   Copy Number Variation (CNV). 1000 predictor features for each subtype is included to satisfy package size limitations.}
#'   \item{Cell_line_Index_Dream}{List of Cell Line names of each genomic charcterization}
#'   \item{finalY_train_Dream}{Drug Sensitivity of training samples (35) for 31 drugs provided for 
#'   NCI-Dream Drug Sensitivity Prediction Challenge }
#'   \item{finalY_train_cell_Dream}{Cell line names of the training samples}
#'   \item{finalY_test_Dream}{Drug Sensitivity of testing samples (18) for 31 drugs provided for 
#'   NCI-Dream Drug Sensitivity Prediction Challenge Dataset}
#'   \item{finalY_test_cell_Dream}{Cell line names of the testing samples}
#' }
#' @source \url{https://www.synapse.org/#!Synapse:syn2785778/wiki/70252}
#' @references
#' Costello, James C., et al. "A community effort to assess and improve drug sensitivity prediction algorithms." Nature biotechnology 32.12 (2014): 1202-1212.
#' 
"Dream_Dataset"
