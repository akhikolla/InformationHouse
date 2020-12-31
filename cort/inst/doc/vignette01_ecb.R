## ----setup, include = FALSE---------------------------------------------------
library(cort)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 7
)

## ----fig.cap="Pairs-plot of original peusdo-observation from the data"--------
set.seed(1)
data("LifeCycleSavings")
pseudo_data <- (apply(LifeCycleSavings,2,rank,ties.method="max")/(nrow(LifeCycleSavings)+1))

pairs(pseudo_data,lower.panel=NULL)

## -----------------------------------------------------------------------------
(cop <- cbCopula(x = pseudo_data,m = 5,pseudo = TRUE))

## ---- fig.cap = "Pairs-plot of original peusdo-observation from the data (red) with simulated pseudo_observation (black)"----
simu <- rCopula(n = 1000,copula = cop)

pairs(rbind(simu,pseudo_data),
      col=c(rep("black",nrow(simu)),rep("red",nrow(pseudo_data))),
      gap=0,
      lower.panel=NULL,cex=0.5)

