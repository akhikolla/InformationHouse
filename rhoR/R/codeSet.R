#' @title codeSet
#' 
#' @description A codeSet is a Nx2 binary matrix in which the first column corresponds to the first rater and the second column corresponds to the second rater.
#' 
#' @format The codeSet is an object of class matrix with n rows and two columns. 
#' 
#' @examples 
#' #An example codeSet
#' firstRater = c(1,1,1,1,rep(0,36))
#' secondRater = c(1,1,1,0,1,1,rep(0,34))
#' exampleSet = cbind(firstRater,secondRater)
#'
#' #This set is included in the package under the variable name "codeSet".
"codeSet"
