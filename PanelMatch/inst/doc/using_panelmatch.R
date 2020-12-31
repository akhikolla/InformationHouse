## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(PanelMatch)
DisplayTreatment(unit.id = "wbcode2",
                 time.id = "year", legend.position = "none",
                 xlab = "year", ylab = "Country Code",
                 treatment = "dem", data = dem)

## -----------------------------------------------------------------------------
DisplayTreatment(unit.id = "wbcode2",
                 time.id = "year", legend.position = "none",
                 xlab = "year", ylab = "Country Code",
                 treatment = "dem", data = dem, 
                 hide.x.axis.label = TRUE, hide.y.axis.label = TRUE) # axis label options

## -----------------------------------------------------------------------------
DisplayTreatment(unit.id = "wbcode2",
                 time.id = "year", legend.position = "none",
                 xlab = "year", ylab = "Country Code",
                 treatment = "dem", data = dem, 
                 hide.x.axis.label = TRUE, hide.y.axis.label = TRUE, 
                 dense.plot = TRUE) # setting dense.plot to TRUE

## -----------------------------------------------------------------------------
PM.results.none <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2", 
                         treatment = "dem", refinement.method = "none", 
                         data = dem, match.missing = TRUE, 
                         size.match = 5, qoi = "att", outcome.var = "y",
                         lead = 0:4, forbid.treatment.reversal = FALSE, 
                         use.diagonal.variance.matrix = TRUE)

## -----------------------------------------------------------------------------
PM.results <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2", 
                         treatment = "dem", refinement.method = "mahalanobis", # use Mahalanobis distance 
                         data = dem, match.missing = TRUE, 
                         covs.formula = ~ tradewb, 
                         size.match = 5, qoi = "att" , outcome.var = "y",
                         lead = 0:4, forbid.treatment.reversal = FALSE, 
                         use.diagonal.variance.matrix = TRUE)

## -----------------------------------------------------------------------------
PM.results <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2", 
                         treatment = "dem", refinement.method = "mahalanobis", 
                         data = dem, match.missing = TRUE, 
                         covs.formula = ~ I(lag(tradewb, 1:4)) + I(lag(y, 1:4)), # lags
                         size.match = 5, qoi = "att", outcome.var = "y",
                         lead = 0:4, forbid.treatment.reversal = FALSE, 
                         use.diagonal.variance.matrix = TRUE)

## -----------------------------------------------------------------------------
PM.results1 <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2", 
                         treatment = "dem", refinement.method = "mahalanobis", 
                         data = dem, match.missing = FALSE, listwise.delete = TRUE, # listwise deletion used 
                         covs.formula = ~ I(lag(tradewb, 1:4)) + I(lag(y, 1:4)), 
                         size.match = 5, qoi = "att", outcome.var = "y",
                         lead = 0:4, forbid.treatment.reversal = FALSE, 
                         use.diagonal.variance.matrix = TRUE)

## -----------------------------------------------------------------------------
PM.results2 <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2", 
                         treatment = "dem", refinement.method = "ps.weight", 
                         data = dem, match.missing = FALSE, listwise.delete = TRUE, 
                         covs.formula = ~ I(lag(tradewb, 1:4)) + I(lag(y, 1:4)), 
                         size.match = 5, qoi = "att", outcome.var = "y",
                         lead = 0:4, forbid.treatment.reversal = FALSE, 
                         use.diagonal.variance.matrix = TRUE)

## -----------------------------------------------------------------------------
# extract the first matched set
mset <- PM.results.none$att[1]

DisplayTreatment(unit.id = "wbcode2",
                 time.id = "year", legend.position = "none",
                 xlab = "year", ylab = "Country Code",
                 treatment = "dem", data = dem,
                 matched.set = mset, # this way we highlight the particular set
                 show.set.only = TRUE)


## -----------------------------------------------------------------------------
summary(PM.results.none$att)

plot(PM.results.none$att)

plot(PM.results.none$att, include.empty.sets = TRUE) # The red tiny bar that would otherwise indicate empty sets is now part of the grey bar

## -----------------------------------------------------------------------------
# get covariate balance for sets that are unrefined

get_covariate_balance(PM.results.none$att,
                      data = dem,
                      covariates = c("tradewb", "y"),
                      plot = FALSE)

# compare with sets that have had various refinements applied

get_covariate_balance(PM.results$att,
                      data = dem,
                      covariates = c("tradewb", "y"),
                      plot = FALSE)

get_covariate_balance(PM.results1$att,
                      data = dem,
                      covariates = c("tradewb", "y"), 
                      plot = FALSE)

get_covariate_balance(PM.results2$att,
                      data = dem,
                      covariates = c("tradewb", "y"), 
                      plot = FALSE)



## -----------------------------------------------------------------------------
get_covariate_balance(PM.results1$att,
                      data = dem,
                      covariates = c("tradewb", "y"), 
                      plot = TRUE, # visualize by setting plot to TRUE
                      ylim = c(-.2, .2))

## -----------------------------------------------------------------------------
balance_scatter(non_refined_set = PM.results.none$att,
               refined_list = list(PM.results$att, PM.results2$att),
               data = dem,
               covariates = c("y", "tradewb"))

## -----------------------------------------------------------------------------
PE.results <- PanelEstimate(sets = PM.results1, data = dem)

names(PE.results)

# View the point estimates 

PE.results[["estimates"]]


## -----------------------------------------------------------------------------
summary(PE.results)

plot(PE.results)

