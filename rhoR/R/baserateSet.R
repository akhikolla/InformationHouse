###
#' @title Calculate Baserate (Set)
#' 
#' @description
#' This function will calculate the baserate of the first rater, second rater, and the average of both the raters. Called by \code{\link{baserate}}.
#' 
#' @param set The \code{\link{codeSet}} for which the baserate is calculated
#' 
#' @seealso \code{\link{baserate}} and \code{\link{baserateCT}}
#' 
#' @export
#' @return A list of the format:\describe{
#' \item{firstBaserate}{The percentage that a positive code, or a 1, appears in the first rater}
#' \item{secondBaserate}{The percentage that a positive code, or a 1, appears in the second rater}
#' \item{averageBaserate}{The average percentage that a positive code, or a 1, appears in either of the two raters}
#' }
#' 
###
baserateSet = function(set){
  firstBaserate = length(which(set[,1] == 1))/nrow(set)
  secondBaserate = length(which(set[,2] == 1))/nrow(set)
  averageBaserate = mean(c(firstBaserate,secondBaserate))
  return(list(firstBaserate = firstBaserate, secondBaserate = secondBaserate, averageBaserate = averageBaserate))
}