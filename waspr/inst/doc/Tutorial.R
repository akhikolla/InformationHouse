## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----message=FALSE------------------------------------------------------------
library(waspr)

## -----------------------------------------------------------------------------

out <- wasp(pois_logistic,
            iter = 10,
            acc = 0.001,
            par.names = c("beta_s", "alpha_l", "beta_l",
                         "baseline_sigma", "baseline_mu",
                         "correlation", "sigma_s", "sigma_l"))


## -----------------------------------------------------------------------------
summary(out)

