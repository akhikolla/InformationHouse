#' @title contingencyTable
#' 
#' @description A contingency Table is a 2x2 matrix that contains the counts of all combinations of positive and negative ratings made by two raters.
#' 
#' @format The contingency Table is an object of class matrix with two rows and two columns.  The ordering of the combination vector input to the matrix is as follows: c(Rater1Positive & Rater2Positive, Rater1Negative & Rater2Positive, Rater1Positive & Rater2Negative, Rater1Negative & Rater2Negative).
#' 
#' @examples 
#' #An example contingencyTable
#' ct = matrix(c(3,2,1,34), nrow = 2, ncol = 2)
#' 
#' #This contingencyTable is included in the package under the variable name "contingencyTable".
#' 
"contingencyTable"