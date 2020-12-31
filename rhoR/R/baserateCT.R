###
#' @title Calculate Baserate (CT)
#' 
#' @description
#' This function calculates the baserate of the first rater, second rater, and the average of both the raters. Called by \code{\link{baserate}}.
#' 
#' @export
#' 
#' @param CT The \code{\link{contingencyTable}} for which the baserate is calculated
#' 
#' @seealso \code{\link{baserate}} and \code{\link{baserateSet}}
#' 
#' @return A list of the format:\describe{
#' \item{firstBaserate}{The percentage of the data for which a positive code, or a 1, appears in the first rater}
#' \item{secondBaserate}{The percentage of the data for which a positive code, or a 1, appears in the second rater}
#' \item{averageBaserate}{The average of the firstBaserate and secondBaserate.}
#' }
#' 
###

baserateCT = function(CT){
  if(any(CT < 0)){stop("Values in Contingency Table must be positive")}
  firstBaserate = sum(CT[1,])/sum(CT)
  secondBaserate = sum(CT[,1])/sum(CT)
  averageBaserate = mean(c(firstBaserate,secondBaserate))
  return(list(firstBaserate = firstBaserate, secondBaserate = secondBaserate, averageBaserate = averageBaserate))
}