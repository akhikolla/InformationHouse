## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
library(surveysd)

set.seed(1234)
eusilc <- demo.eusilc(prettyNames = TRUE)
dat_boot <- draw.bootstrap(eusilc, REP = 10, hid = "hid", weights = "pWeight",
                           strata = "region", period = "year")
dat_boot_calib <- recalib(dat_boot, conP.var = "gender", conH.var = "region",
                          epsP = 1e-2, epsH = 2.5e-2, verbose = FALSE)
dat_boot_calib[, onePerson := nrow(.SD) == 1, by = .(year, hid)]

## print part of the dataset
dat_boot_calib[1:5, .(year, povertyRisk, eqIncome, onePerson, pWeight, w1, w2, w3, w4, w5)]

## -----------------------------------------------------------------------------
povertyRate <- calc.stError(dat_boot_calib, var = "povertyRisk", fun = weightedRatio)
totalIncome <- calc.stError(dat_boot_calib, var = "eqIncome", fun = weightedSum)

## -----------------------------------------------------------------------------
povertyRate$Estimates
totalIncome$Estimates

## -----------------------------------------------------------------------------
## define custom estimator
myWeightedSum <- function(x, w) {
  sum(x*w)
}

## check if results are equal to the one using `suveysd::weightedSum()`
totalIncome2 <- calc.stError(dat_boot_calib, var = "eqIncome", fun = myWeightedSum)
all.equal(totalIncome$Estimates, totalIncome2$Estimates)

## -----------------------------------------------------------------------------
multipleRates <- calc.stError(dat_boot_calib, var = c("povertyRisk", "onePerson"), fun = weightedRatio)
multipleRates$Estimates

## -----------------------------------------------------------------------------
dat2 <- subset(dat_boot_calib, year == 2010)
for (att  in c("period", "weights", "b.rep"))
  attr(dat2, att) <- attr(dat_boot_calib, att)

## -----------------------------------------------------------------------------
povertyRates <- calc.stError(dat2, var = "povertyRisk", fun = weightedRatio, group = "region")
povertyRates$Estimates

## -----------------------------------------------------------------------------
povertyRates <- calc.stError(dat2, var = "povertyRisk", fun = weightedRatio, 
                             group = c("gender", "region"))
povertyRates$Estimates


## -----------------------------------------------------------------------------
povertyRates <- calc.stError(dat2, var = "povertyRisk", fun = weightedRatio, 
                             group = list(c("gender", "region")))
povertyRates$Estimates

## -----------------------------------------------------------------------------
povertyRates <- calc.stError(dat2, var = "povertyRisk", fun = weightedRatio, 
                             group = list("gender", "region", c("gender", "region")))
povertyRates$Estimates

