plot.outlierProbs <- function(x,...) {
  dotchart(as.numeric(x), xlim = c(0,1),xlab="Outlier Probability",...)
  abline(v=0.0,lty=3)
  abline(v=1.0,lty=3)
}