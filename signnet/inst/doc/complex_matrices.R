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

## ----projection_example-------------------------------------------------------
# construct network
el <- matrix(c(1,"a",1,"b",1,"c",2,"a",2,"b"),ncol = 2,byrow = TRUE)
g <- graph_from_edgelist(el,directed = FALSE)
E(g)$sign <- c(1,1,-1,1,-1)
V(g)$type <- c(FALSE,TRUE,TRUE,TRUE,FALSE)

# vertex duplication
gu <- as_unsigned_2mode(g,primary = TRUE)

# project and binarize
pu <- bipartite_projection(gu,which = "true")
pu <- delete_edge_attr(pu,"weight")

# vertex contraction
ps <- as_signed_proj(pu)
igraph::as_data_frame(ps,"edges")

## ----pna_adj------------------------------------------------------------------
as_adj(ps,type = "both", attr = "type", sparse = FALSE)

## ----complex_adj,eval=FALSE---------------------------------------------------
#  as_adj_complex(ps,attr = "type")

## ----sneaky_show,echo=FALSE---------------------------------------------------
structure(c(0+0i, 0.5-0.5i, 0-1i, 0.5+0.5i, 0+0i, 0-1i, 0+1i, 
0+1i, 0+0i), .Dim = c(3L, 3L), .Dimnames = list(c("a", "b", "c"
), c("a", "b", "c")))

## ----complex_lapl, eval=FALSE-------------------------------------------------
#  laplacian_matrix_complex(ps, attr = "type")

## ----sneaky_show1,echo=FALSE--------------------------------------------------
structure(c(1.70710678118655+0i, -0.5+0.5i, 0+1i, -0.5-0.5i, 
1.70710678118655+0i, 0+1i, 0-1i, 0-1i, 2+0i), .Dim = c(3L, 3L
), .Dimnames = list(c("a", "b", "c"), c("a", "b", "c")))

## ----ambi_net-----------------------------------------------------------------
g <- graph.full(5)
E(g)$type <- c(rep("P",3),rep("N",3),rep("A",4))

count_complex_triangles(g,attr = "type")

