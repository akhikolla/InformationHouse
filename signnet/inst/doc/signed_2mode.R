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

## ----simple s2mode------------------------------------------------------------
el <- matrix(c(1,"a",1,"b",1,"c",2,"a",2,"b"),ncol = 2,byrow = TRUE)
g <- graph_from_edgelist(el,directed = FALSE)
E(g)$sign <- c(1,1,-1,1,-1)
V(g)$type <- c(FALSE,TRUE,TRUE,TRUE,FALSE)

## ----graph, echo=FALSE,out.width = "50%",fig.align='center'-------------------
knitr::include_graphics("small_signed2mode.png")

## ----matmul-------------------------------------------------------------------
A <- as_incidence_signed(g)
R <- A%*%t(A)
C <- t(A)%*%A
R
C

## ----duplicate----------------------------------------------------------------
gu <- as_unsigned_2mode(g,primary = TRUE)
gu

## ----binarize-----------------------------------------------------------------
pu <- bipartite_projection(gu,which = "true")
pu <- delete_edge_attr(pu,"weight")
pu

## ----contract-----------------------------------------------------------------
ps <- as_signed_proj(pu)
as_data_frame(ps,"edges")

