## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message=FALSE,
  warning=FALSE
)

## ----setup--------------------------------------------------------------------
library(igraph)
library(signnet)

## ----deg_same-----------------------------------------------------------------
A <- matrix(c(0,  1,  0,  1,  0,  0,  0, -1, -1,  0,  
               1,  0,  1, -1,  1, -1, -1,  0,  0,  0,  
               0,  1,  0,  1, -1,  0,  0,  0, -1,  0,  
               1, -1,  1,  0,  1, -1, -1,  0,  0,  0,  
               0,  1, -1,  1,  0,  1,  0, -1,  0, -1,  
               0, -1,  0, -1,  1,  0,  1,  0,  1, -1,  
               0, -1,  0, -1,  0,  1,  0,  1, -1,  1,  
              -1,  0,  0,  0, -1,  0,  1,  0,  1,  0,  
              -1,  0, -1,  0,  0,  1, -1,  1,  0,  1,  
               0,  0,  0,  0, -1, -1,  1,  0,  1,  0),10,10)

g <- graph_from_adjacency_matrix(A,"undirected",weighted = "sign")

degree_signed(g,type="ratio")
eigen_centrality_signed(g)
pn_index(g)


## ----cor----------------------------------------------------------------------
cor(eigen_centrality_signed(g),pn_index(g),method = "kendall")

## ----non_dominant,error=TRUE--------------------------------------------------
A <- matrix(c( 0,  1,  1, -1,  0,  0, -1,  0,  0, 
               1,  0,  1,  0, -1,  0,  0, -1,  0, 
               1,  1,  0,  0,  0, -1,  0,  0, -1, 
              -1,  0,  0,  0,  1,  1, -1,  0,  0, 
               0, -1,  0,  1,  0,  1,  0, -1,  0, 
               0,  0, -1,  1,  1,  0,  0,  0, -1, 
              -1,  0,  0, -1,  0,  0,  0,  1,  1, 
               0, -1,  0,  0, -1,  0,  1,  0,  1, 
               0,  0, -1,  0,  0, -1,  1,  1, 0), 9, 9)

g <- igraph::graph_from_adjacency_matrix(A,"undirected",weighted = "sign")
eigen_centrality_signed(g)


