#' parse_input.R
#'
#' Internal function to accept the input as unsorted Z-statistics and either a matrix
#' or vector of correlations, and return the t_vec and pairwise_cors. Also do limited
#' error checking.
#'
#' @param test_stats All test statistics in the set
#' @param cor_mat The correlation matrix of the test statistics
#' @param pairwise_cors The vector of all pairwise correlations
#'
#' @return A list with the elements:
#' \item{t_vec}{The sorted magnitudes of test statistics}
#' \item{pairwise_cors}{The vector of all pairwise correlations}
#'
#' @keywords internal
#' @export
#' @examples
#' parse_input(test_stats=rnorm(5), pairwise_cors=rep(0.3,10))

parse_input <- function(test_stats, cor_mat=NULL, pairwise_cors=NULL) {

  # Limit the size of sets to 900 factors
  if (length(test_stats) > 2000) {stop('You have too many factors, please restrict to 2000')}

	# Ensure that the thresholds are sorted in descending order, largest first.
	t_vec <- sort(abs(test_stats), decreasing=TRUE)
	d <- length(t_vec)

	# Sometimes test stats are too big for R's precision
	too_big <- which(t_vec > 8.2)
	if (length(too_big) > 0) {t_vec[too_big] <- 8.2}

	# Did they specify correlations?
	if (is.null(cor_mat) & is.null(pairwise_cors)) {
		stop("You must specify either cor_mat or pairwise_cors!")
	}

	# cor_mat specification gets priority
	if (!is.null(cor_mat)) {
		if(!isSymmetric(cor_mat)) {
			stop("You did not specify a symmetric correlation matrix")
		}

	  # Put the correlation matrix into pairwise_cor
	  pairwise_cors <- cor_mat[upper.tri(cor_mat)]
	}

	# Correct number of pairwise correlations?
	if (length(pairwise_cors) != d*(d-1)/2) {
		stop("Your pairwise correlation matrix/vector is of the wrong size!")
	}

	return( list(t_vec=t_vec, pairwise_cors=pairwise_cors))
}
