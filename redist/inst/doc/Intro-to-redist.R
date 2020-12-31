## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(redist)
library(igraph)
library(spdep)
library(coda)
set.seed(1)

## ---- eval=F------------------------------------------------------------------
#  install.packages('redist')

## ---- eval=F------------------------------------------------------------------
#  devtools::install_github(repo = 'kosukeimai/redist', ref = 'master')

## -----------------------------------------------------------------------------
library(sf)
library(sp)

## -----------------------------------------------------------------------------
library(spdep)

## -----------------------------------------------------------------------------
library(igraph)
library(ggplot2)

## -----------------------------------------------------------------------------
library(dplyr)
library(magrittr)

## ----data---------------------------------------------------------------------
data("algdat.p10")
data("algdat.p20")
data("algdat.pfull")

## ----objects------------------------------------------------------------------
names(algdat.p10)
names(algdat.p20)
names(algdat.pfull)

## -----------------------------------------------------------------------------
data("fl25")

## -----------------------------------------------------------------------------
head(fl25)

## -----------------------------------------------------------------------------
data("fl25")

## ---- fig.width=5-------------------------------------------------------------
fl25 %>% ggplot(aes(fill = pop)) + 
  geom_sf()

## ---- fig.width=5-------------------------------------------------------------
fl25$id <- 1:25
fl25[c(18,19,23:25),] %>% ggplot() + 
  geom_sf() +
  geom_sf_label(aes(label = id))

## -----------------------------------------------------------------------------
adjlist <- poly2nb(pl = fl25, queen = FALSE)

## -----------------------------------------------------------------------------
adjlist[[25]]

## ---- fig.width=5-------------------------------------------------------------
plot(graph_from_adj_list(adjlist, mode = 'total'))

## -----------------------------------------------------------------------------
for(i in 1:25){
  adjlist[[i]] <- adjlist[[i]]-1
}

## -----------------------------------------------------------------------------
adjlist[[25]]

## -----------------------------------------------------------------------------
alg_mcmc <- redist.mcmc(adjobj = algdat.pfull$adjlist,
                        popvec = algdat.pfull$precinct.data$pop,
                        ndists = 3,
                        nsims = 10000,
                        savename = "redist.mcmc")

## -----------------------------------------------------------------------------
initcds <- algdat.pfull$cdmat[,1]

## -----------------------------------------------------------------------------
alg_mcmc <- redist.mcmc(adjobj = algdat.pfull$adjlist,
                        popvec = algdat.pfull$precinct.data$pop,
                        initcds = initcds,
                        nsims = 10000,
                        savename = "redist.mcmc")

## -----------------------------------------------------------------------------
class(alg_mcmc)
names(alg_mcmc)

## -----------------------------------------------------------------------------
alg_mcmc$partitions[,1]

## -----------------------------------------------------------------------------
alg_mcmc$distance_parity[1]

## ---- eval=F------------------------------------------------------------------
#  library(Rmpi)
#  redist.mcmc.mpi(adjobj = algdat.pfull$adjlist,
#                  popvec = algdat.pfull$precinct.data$pop,
#                  nsims = 10000,
#                  ndists = 3,
#                  savename = "redist.mcmc.mpi")

## -----------------------------------------------------------------------------
RNGkind(kind = "L'Ecuyer-CMRG")
set.seed(1)
nchains <- 4
nsims <- 10000

## -----------------------------------------------------------------------------
mcmc_chains <- lapply(1:nchains, function(x){
          redist.mcmc(adjobj = algdat.pfull$adjlist, 
                      popvec = algdat.pfull$precinct.data$pop, 
                      nsims = nsims,
                      ndists = 3)
})

## ---- eval=F------------------------------------------------------------------
#  mcmc_chains <- parallel::mclapply(1:nchains, function(x){
#            redist.mcmc(adjobj = algdat.pfull$adjlist,
#                        popvec = algdat.pfull$precinct.data$pop,
#                        nsims = nsims,
#                        ndists = 3)
#  }, mc.set.seed = 1, mc.cores = parallel::detectCores())

## -----------------------------------------------------------------------------
rsg <- redist.rsg(adj.list = algdat.pfull$adjlist,
                  population = algdat.pfull$precinct.data$pop,
                  ndists = 3,
                  thresh = 0.05, 
                  maxiter = 5000)

## -----------------------------------------------------------------------------
rsg$district_membership

## -----------------------------------------------------------------------------
rsg$district_list

## -----------------------------------------------------------------------------
rsg$district_pop

## -----------------------------------------------------------------------------
seg <- redist.segcalc(algout = alg_mcmc, 
                      grouppop = algdat.pfull$precinct.data$blackpop,
                      fullpop = algdat.pfull$precinct.data$pop)

## -----------------------------------------------------------------------------
redist.diagplot(seg, plot = "trace")

## -----------------------------------------------------------------------------
redist.diagplot(seg, plot = "autocorr")

## -----------------------------------------------------------------------------
redist.diagplot(seg, plot = "densplot")

## -----------------------------------------------------------------------------
redist.diagplot(seg, plot = "mean")

## -----------------------------------------------------------------------------
seg_chains <- lapply(1:nchains, 
                     function(i){redist.segcalc(algout = mcmc_chains[[i]], 
                                  grouppop = algdat.pfull$precinct.data$blackpop,
                                  fullpop = algdat.pfull$precinct.data$pop)})

redist.diagplot(seg_chains, plot = 'gelmanrubin')

