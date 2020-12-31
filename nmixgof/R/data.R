#' Northern shoveler data
#'
#' Repeated count data of Northern shoveler with covariates, 
#' formatted for use with the \link[unmarked]{unmarked} package.
#'
#' @docType data
#' 
#' @format A list with three elements
#' \describe{
#'   \item{y}{A matrix with Northern shoveler counts}
#'   \item{site}{A data frame with site specific covariates}
#'   \item{obs}{A list containing observation specific covariates}
#' }
#' 
#' @references Knape et al. (2018) Methods in Ecology and Evolution in press.
#' (\href{https://www.biorxiv.org/content/early/2017/09/27/194340}{BioRxiv})
#' 
#' @examples
#' library(unmarked)
#' umf = unmarkedFramePCount(y = shoveler$y, obsCovs = shoveler$obs, siteCovs = shoveler$site)
"shoveler"