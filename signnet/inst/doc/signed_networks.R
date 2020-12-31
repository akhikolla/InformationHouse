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

## ----example------------------------------------------------------------------
g <- graph.full(5,directed = FALSE,loops = FALSE)
E(g)$sign <- 1
g

## ----signed_adj---------------------------------------------------------------
data("tribes")
as_adj_signed(tribes)[1:5,1:5]

## ----signed_lap---------------------------------------------------------------
laplacian_matrix_signed(tribes)[1:5,1:5]

