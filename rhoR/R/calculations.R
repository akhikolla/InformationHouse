rating_set_recall <- function (x) {
  if (!is(x, "contingency.table")) {
    x <- as.contingency.table(x)
  }
  
  x$agreement$true_positive / (x$agreement$true_positive + x$agreement$false_negative)
}

rating_set_precision <- function (x) {
  if (!is(x, "contingency.table")) {
    x <- as.contingency.table(x)
  }
  
  x$agreement$true_positive / (x$agreement$true_positive + x$agreement$false_positive)
}