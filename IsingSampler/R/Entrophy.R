IsingEntrophy <- function(
  graph,
  thresholds,
  beta = 1,
  conditional = numeric(0), # Indices of nodes to condition on
  marginalize = numeric(0),
  base = 2,
  responses = c(0L, 1L)
  ){
  stopifnot(isSymmetric(graph))  
  stopifnot(length(responses)==2)
  if (any(diag(graph)!=0))
  {
    diag(graph) <- 0
    warning("Diagonal set to 0")
  }

  Lik <- IsingLikelihood(graph, thresholds, beta, responses)

  varNames <- names(Lik)[-1]

  if (any(marginalize %in% conditional)){
    stop("can not marginalize over nodes to condition on")
  }
  
  if (length(marginalize) > 0){
    Lik <- Lik %>% group_by_(.dots = varNames[-marginalize]) %>%
      dplyr::summarize_(Probability = ~sum(Probability))
  }
  
  if (length(conditional) > 0){
    Lik <- Lik %>% group_by_(.dots = varNames[conditional])
  } else {
    Lik <- Lik %>% ungroup()
  }
  
  condLik <- Lik %>% 
    dplyr::summarize_(
      P = ~sum(Probability),
      Entrophy = ~-sum(Probability/sum(Probability) * log(Probability/sum(Probability), base) )
      )
  
  Ent <- sum(condLik$P * condLik$Entrophy)
  
  return(Ent)
}


NodeInformation <- function(
  graph,
  thresholds,
  beta = 1,
  base = 2,
  responses = c(0L, 1L)
){
  stopifnot(isSymmetric(graph))  
  stopifnot(length(responses)==2)
  if (any(diag(graph)!=0))
  {
    diag(graph) <- 0
    warning("Diagonal set to 0")
  }
  
  # Shannon information of node to whole graph
  
  sapply(seq_len(ncol(graph)),
         function(Node){
           graphEntrophy <- IsingEntrophy(graph, thresholds, beta, responses=responses,base=base, marginalize = Node)
           nodeEntrophy <- IsingEntrophy(graph, thresholds, beta, responses=responses,base=base, conditional = Node)
           return(graphEntrophy - nodeEntrophy)
         })
}
