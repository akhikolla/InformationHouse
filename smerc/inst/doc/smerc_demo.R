## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
library(smerc) # load package
data(nydf) # load data
str(nydf) # look at structure

## -----------------------------------------------------------------------------
data(nypoly) # load nypoly data
library(sp) # load sp package for plotting
plot(nypoly) # plot SpatialPolygonsDataFrame

## -----------------------------------------------------------------------------
coords = nydf[,c("x", "y")] # extract coordinates
cases = nydf$cases # extract cases
pop = nydf$population # extract population
scan_out = scan.test(coords, cases, pop, nsim = 99) # perform scan test

## -----------------------------------------------------------------------------
class(scan_out)

## -----------------------------------------------------------------------------
scan_out # print scan.test results

## -----------------------------------------------------------------------------
summary(scan_out) # summarize scan.test results

## -----------------------------------------------------------------------------
plot(scan_out) # basic plot of scan.test results

## -----------------------------------------------------------------------------
plot(nypoly, col = color.clusters(scan_out)) #nicer plot of scan.test results

## -----------------------------------------------------------------------------
bn_out = bn.test(coords = coords, cases = cases, pop = pop, cstar = 20,
                 alpha = 0.01) # perform besag-newell test
bn_out # print results
summary(bn_out) # summarize results
plot(bn_out) # plot results

## -----------------------------------------------------------------------------
# perform CEPP test
cepp_out = cepp.test(coords = coords, cases = cases, pop = pop,
                     nstar = 5000, nsim = 99, alpha = 0.1)
cepp_out # print results
summary(cepp_out) # summarize results
plot(cepp_out) # plot results

## -----------------------------------------------------------------------------
w = dweights(coords, kappa = 1) # construct weights matrix
tango_out = tango.test(cases, pop, w, nsim = 49) # perform tango's test
tango_out # print results
plot(tango_out) # plot results

## -----------------------------------------------------------------------------
# obtain zones for elliptical scan method
ezones = elliptic.zones(coords, pop, ubpop = 0.1)
# view structure of ezones
str(ezones)

