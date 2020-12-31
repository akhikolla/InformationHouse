# Wrapper function for:
# - pl: pseudolikelihood
# - uni: univariate logistic regressions
# - bi: bivariate logistic regressions
# - ll: Loglinear model

EstimateIsing <- function(data, responses, beta = 1, method = c('pl', 'uni', 'bi', 'll'),adj = matrix(1, ncol(data), ncol(data)), ...){

  method <- match.arg(method)
  
  if (!identical(adj,matrix(1, ncol(data), ncol(data)))){

    if (!method %in% c("uni","ll")) stop("Adjacency structure only supported if method = 'll' or method = 'uni'")
  }

  data <- as.matrix(data)
 
  switch(method[[1]],
        pl = EstimateIsingPL(data, responses, beta, ...)  ,
        uni = EstimateIsingUni(data, responses, beta, adj=adj,...)  ,
        bi = EstimateIsingBi(data, responses, beta, ...),
        ll = EstimateIsingLL(data, responses, beta, adj=adj,...)
         )
}