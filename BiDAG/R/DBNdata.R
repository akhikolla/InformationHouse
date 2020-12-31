#' A simulated data set from a 2-step dynamic Bayesian network
#'  
#' A synthetic dataset containing 100 observations generated from a random dynamic Bayesian network with 12 continuous dynamic nodes and 3 static discrete nodes.
#' The DBN imcludes observations from 5 time slices.
#'
#' @format A data frame with 100 rows and 63 (3+12*5) columns representing observations of 15 variables: 3 static variables (first 3 columns) which do not change over time and 12 dynamic variables observed in 5 conseecutive time slices.
#'
"DBNdata"
