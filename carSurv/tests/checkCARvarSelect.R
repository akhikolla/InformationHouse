##########################################################
# Compare weighted unbiased variance C version with R code
library(carSurv)

# Is the value NA, if wrong values are passed? 
checkRes <- carVarSelect(carSurvScores=c(NA, 5, Inf, 6:10))
stopifnot(is.na(checkRes))
