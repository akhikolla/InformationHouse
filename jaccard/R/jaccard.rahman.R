#' Compute p-value using an extreme value distribution
#'
#' Rahman et al. (2014) proposes a method to compute a p-value of 
#' a Jaccard/Tanimoto coefficient using an extreme value distribution.
#' Their paper provides the following description:
#' The mean (mu) and s.d. (sigma) of the similarity scores are used to
#' define the z score, z = (Tw - mu)/sigma. For the purpose of calculating
#' the P value, only hits with T > 0 are considered. The P value w
#' is derived from the z score using an extreme value distribution
#' P = 1 - exp(-e-z*pi/sqrt(6) - G'(1)), where the Euler=Mascheroni constant G'(1)=0.577215665.
#'
#' @param j a numeric vector of observed Jaccard coefficients (uncentered)
#' @return \code{jaccard.rahman} returns a numeric vector of p-values
#'
#' @references Rahman, Cuesta, Furnham, Holliday, and Thornton (2014) EC-BLAST: a tool to automatically search and compare enzyme reactions. Nature Methods, 11(2) \url{http://www.nature.com/nmeth/journal/v11/n2/full/nmeth.2803.html}
#'
#' @export jaccard.rahman
jaccard.rahman <- function(j) {
	mu = mean(j)
	sd = sd(j)
	z = (j - mu) / sd
	p = 1 - exp(-exp(-z*pi/sqrt(6)-0.577215665))
	return(p)
}
