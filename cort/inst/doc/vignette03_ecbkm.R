## ----setup, include = FALSE---------------------------------------------------
library(cort)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 7
)

## ----fig.cap="Pairs-plot of original peusdo-observations"---------------------
set.seed(1)
data("LifeCycleSavings")
dataset <- (apply(LifeCycleSavings,2,rank,ties.method="max")/(nrow(LifeCycleSavings)+1))
pairs(dataset,col="2",lower.panel=NULL)


## -----------------------------------------------------------------------------
  known_margins <- c(2,3)
  known_copula <- cbCopula(x = dataset[,known_margins],m = 25,pseudo = TRUE)

## -----------------------------------------------------------------------------
  (cop <- cbkmCopula(x = dataset,m = 5,pseudo = TRUE,margins_numbers = known_margins,known_cop = known_copula))

## ---- fig.cap="Pairs-plot of the original data (red) and simulated data from a good model (black)"----
  simu <- rCopula(1000,cop)
  pairs(rbind(simu,dataset),col=c(rep("black",nrow(simu)),rep("red",nrow(dataset))),gap=0,lower.panel = NULL,cex=0.5)

## ---- fig.cap="Pairs-plot of the original data (red) and simulated data from a wrong model (black)"----
  bad_known_copula <- cbCopula(x = dataset[,known_margins],m = 2,pseudo = TRUE)
  cop <- cbkmCopula(x = dataset,m = 5,pseudo = TRUE,margins_numbers = known_margins,known_cop = bad_known_copula)

  simu <- rCopula(1000,cop)
  pairs(rbind(simu,dataset),col=c(rep("black",nrow(simu)),rep("red",nrow(dataset))),gap=0,lower.panel = NULL,cex=0.5)

## ---- fig.cap="Pairs-plot of the original data with independance"-------------

  true_copula1 <- known_copula
  true_copula2 <- bad_known_copula

  dataset <- cbind(rCopula(1000,true_copula1),rCopula(1000,true_copula2))
  colnames(dataset) <- c("u","v","w","x")
  pairs(dataset,lower.panel=NULL,cex=0.5)


## ---- fig.cap="Pairs-plot of the original data (red) and simulated data (black) -- Independance case"----
  cop <- cbkmCopula(x = dataset,m = 5,pseudo = TRUE,margins_numbers = c(1,2),known_cop = true_copula1)
  simu <- rCopula(500,cop)
  pairs(rbind(simu,dataset),col=c(rep("black",nrow(simu)),rep("red",nrow(dataset))),gap=0,lower.panel = NULL,cex=0.5)

