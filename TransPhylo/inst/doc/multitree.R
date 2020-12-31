## ----echo=F--------------------------------------------------------------
knitr::opts_chunk$set(fig.width=6, fig.height=4) 

## ---- library, message=FALSE---------------------------------------------
library(ape)
library(TransPhylo)

## ----simulate, results='hide'--------------------------------------------
neg <- 100/365
off.r <- 1.5
w.shape <- 10
w.scale <- 0.1
ws.shape <- w.shape
ws.scale <- w.scale
pi <- 0.8

set.seed(1234)
simu1 <- simulateOutbreak(neg=neg, off.r=off.r, pi=pi, w.shape=w.shape,
                         w.scale=w.scale, dateStartOutbreak=2000,dateT=2005)
simu2 <- simulateOutbreak(neg=neg, off.r=off.r, pi=pi, w.shape=w.shape,
                         w.scale=w.scale, dateStartOutbreak=2000,dateT=2005)

## ------------------------------------------------------------------------
plot(simu1)
plot(simu2)

## ----results='hide'------------------------------------------------------
ptree1 <- extractPTree(simu1)
ptree2 <- extractPTree(simu2)

iters <- 2e3; thin <- 10
record_tp1 <- inferTTree(ptree1, w.shape, w.scale, ws.shape, ws.scale,
                         mcmcIterations = iters, thinning = thin, dateT = 2005)
record_tp2 <- inferTTree(ptree2, w.shape, w.scale, ws.shape, ws.scale,
                         mcmcIterations = iters, thinning = thin, dateT = 2005)
record_tpj <- infer_multittree_share_param(list(ptree1,ptree2), w.shape, w.scale, ws.shape, ws.scale, 
                                    mcmcIterations = iters, thinning = thin, dateT = 2005,
                                    share = c("neg", "off.r", "off.p", "pi"))

## ------------------------------------------------------------------------
begin <- (iters * 0.5) / thin + 1
end <- iters / thin
record_tp1 <- record_tp1[begin:end]
record_tp2 <- record_tp2[begin:end]
record_tpj[[1]] <- record_tpj[[1]][begin:end]
record_tpj[[2]] <- record_tpj[[2]][begin:end]

## ------------------------------------------------------------------------
get_param_estimates <- function(record, p){
  sapply(record, function(x) x[[p]])
}

df <- data.frame(run = rep(c("tp1","tp2","tp_multitree"), each = length(record_tp1)),
           pi = c(get_param_estimates(record_tp1, "pi"),
                  get_param_estimates(record_tp2, "pi"),
                  get_param_estimates(record_tpj[[1]], "pi")),
           off.r = c(get_param_estimates(record_tp1, "off.r"),
                     get_param_estimates(record_tp2, "off.r"),
                     get_param_estimates(record_tpj[[1]], "off.r")),
           neg = c(get_param_estimates(record_tp1, "neg"),
                   get_param_estimates(record_tp2, "neg"),
                   get_param_estimates(record_tpj[[1]], "neg")))

## ------------------------------------------------------------------------
boxplot(df$pi~df$run,ylab='pi')
lines(c(0,4),rep(0.8,2),col='grey')

## ------------------------------------------------------------------------
boxplot(df$off.r~df$run,ylab='off.r')
lines(c(0,4),rep(1.5,2),col='grey')


## ------------------------------------------------------------------------
boxplot(df$neg~df$run,ylab='Ne*g')
lines(c(0,4),rep(100/365,2),col='grey')

