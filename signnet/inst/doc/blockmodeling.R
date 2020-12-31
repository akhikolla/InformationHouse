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

## ----blockmod_ex--------------------------------------------------------------
g <- sample_islands_signed(10,10,1,20)
clu <- signed_blockmodel(g,k = 10,alpha = 0.5)
table(clu$membership)
clu$criterion

## ----blockmodel_ex_plot,eval=FALSE--------------------------------------------
#  ggblock(g,clu$membership,show_blocks = TRUE)

## ----example, echo=FALSE,out.width = "80%",fig.align='center'-----------------
knitr::include_graphics("blockmodel_example.png")

## ----blockmodel_tribes--------------------------------------------------------
data("tribes")
set.seed(44) #for reproducibility

signed_blockmodel(tribes,k = 3,alpha=0.5,annealing = TRUE)
signed_blockmodel(tribes,k = 3,alpha=0.5,annealing = FALSE)

## ----general_example----------------------------------------------------------
g1 <- g2 <- g3 <- graph.full(5)

V(g1)$name <- as.character(1:5)
V(g2)$name <- as.character(6:10)
V(g3)$name <- as.character(11:15)

g <- Reduce("%u%",list(g1,g2,g3))
E(g)$sign <- 1
E(g)$sign[1:10] <- -1
g <- add.edges(g,c(rbind(1:5,6:10)),attr = list(sign=-1))
g <- add.edges(g,c(rbind(1:5,11:15)),attr = list(sign=-1))
g <- add.edges(g,c(rbind(11:15,6:10)),attr = list(sign=1))

## ----general_blocks-----------------------------------------------------------
set.seed(424) #for reproducibility
blockmat <- matrix(c(1,-1,-1,-1,1,1,-1,1,-1),3,3,byrow = TRUE)
blockmat

general <- signed_blockmodel_general(g,blockmat,alpha = 0.5)
traditional <- signed_blockmodel(g,k = 3,alpha = 0.5,annealing = TRUE)

c(general$criterion,traditional$criterion)

## ----general, echo=FALSE,out.width = "90%",fig.align='center'-----------------
knitr::include_graphics("blockmodel_general.png")

