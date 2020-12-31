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

## ----triangles, echo=FALSE,out.width = "100%"---------------------------------
knitr::include_graphics("balance_triples.png")

## ----example_balanced---------------------------------------------------------
g <- sample_islands_signed(islands.n = 2,islands.size = 10,
                           islands.pin = 0.8,n.inter = 5)

## ----triangle_count-----------------------------------------------------------
count_signed_triangles(g)

## ----triangle_list------------------------------------------------------------
head(signed_triangles(g))

## ----perfect_balance----------------------------------------------------------
balance_score(g, method = "triangles")
balance_score(g, method = "walk")
balance_score(g, method = "frustration")

## ----nonperfect---------------------------------------------------------------
data("tribes")
balance_score(tribes, method = "triangles")
balance_score(tribes, method = "walk")
balance_score(tribes, method = "frustration")

