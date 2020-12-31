## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
library(surveysd)

set.seed(1234)
eusilc <- demo.eusilc(n = 2, prettyNames = TRUE)

eusilc[1:5, .(year, povertyRisk, gender, pWeight)]

## -----------------------------------------------------------------------------
dat_boot <- draw.bootstrap(eusilc, REP = 10, hid = "hid", weights = "pWeight", 
                           strata = "region", period = "year")

## -----------------------------------------------------------------------------
dat_boot_calib <- recalib(dat_boot, conP.var = "gender", conH.var = "region",
                          epsP = 1e-2, epsH = 2.5e-2, verbose = FALSE)
dat_boot_calib[1:5, .(year, povertyRisk, gender, pWeight, w1, w2, w3, w4)]

## -----------------------------------------------------------------------------
err.est <- calc.stError(dat_boot_calib, var = "povertyRisk", fun = weightedRatio, group = "gender")
err.est$Estimates

## -----------------------------------------------------------------------------
group <- list("gender", "region", c("gender", "region"))
err.est <- calc.stError(dat_boot_calib, var = "povertyRisk", fun = weightedRatio, group = group)
head(err.est$Estimates)
## skipping 54 more rows

