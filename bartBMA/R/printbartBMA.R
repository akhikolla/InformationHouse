#write print, summary and predict S3 classes to call for package
#print will take in an object of class "bartBMA"
#' @importFrom utils head
print.bartBMA<-function(x, ...){
  cat("Call:\n")
  print(x$call)
  cat("Training Fitted Values:\n")
  print(c(as.numeric(head(round(x$fitted.values,2)))," ..."))
  cat("Sums of Trees in Occam's Window:\n")
  print(x$sumoftrees)
}