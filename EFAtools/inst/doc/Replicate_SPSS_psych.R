## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.align = "center"
)


if (!requireNamespace("psych", quietly = TRUE)) {
      stop("Package \"psych\" needed for this vignette to work. Please install it.",
      call. = FALSE)
}


## ----setup--------------------------------------------------------------------
library(psych)
library(EFAtools)

## -----------------------------------------------------------------------------
# EFAtools::EFA with type = "psych" without rotation
EFA_psych_paf <- EFA(DOSPERT$cormat, n_factors = 10, N = DOSPERT$N,
                     type = "psych")
# EFAtools::EFA with type = "SPSS" without rotation
EFA_SPSS_paf <- EFA(DOSPERT$cormat, n_factors = 10, N = DOSPERT$N,
                    type = "SPSS")

## -----------------------------------------------------------------------------
# psych::fa without rotation
psych_paf <- psych::fa(DOSPERT$cormat, nfactors = 10, n.obs = DOSPERT$N,
                       fm = "pa", rotate = "none")

## -----------------------------------------------------------------------------
# Compare loadings from psych::fa and EFAtools::EFA with type = "psych"
COMPARE(EFA_psych_paf$unrot_loadings, psych_paf$loadings)

# Compare loadings from SPSS and EFAtools::EFA with type = "SPSS"
COMPARE(EFA_SPSS_paf$unrot_loadings, SPSS$DOSPERT$paf_load)


## -----------------------------------------------------------------------------
# Compare loadings from psych::fa and SPSS
COMPARE(psych_paf$loadings, SPSS$DOSPERT$paf_load)

## -----------------------------------------------------------------------------
## Fit the models

# EFAtools::EFA with type = "psych" with varimax rotation
EFA_psych_var <- EFA(DOSPERT$cormat, n_factors = 10, N = DOSPERT$N,
                     type = "psych", rotation = "varimax")
# EFAtools::EFA with type = "SPSS" with varimax rotation
EFA_SPSS_var <- EFA(DOSPERT$cormat, n_factors = 10, N = DOSPERT$N,
                    type = "SPSS", rotation = "varimax")
# psych::fa with varimax rotation
psych_var <- psych::fa(DOSPERT$cormat, nfactors = 10, n.obs = DOSPERT$N,
                       fm = "pa", rotate = "varimax")

## Check replication of results

# Compare loadings from psych::fa and EFAtools::EFA with type = "psych"
COMPARE(EFA_psych_var$rot_loadings, psych_var$loadings)

# Compare loadings from SPSS and EFAtools::EFA with type = "SPSS"
COMPARE(EFA_SPSS_var$rot_loadings, SPSS$DOSPERT$var_load)

## Compare original results (just to see the difference)

# Compare loadings from psych::fa and SPSS
COMPARE(psych_var$loadings, SPSS$DOSPERT$var_load)


## -----------------------------------------------------------------------------
## Fit the models

# EFAtools::EFA with type = "psych" with promax rotation
EFA_psych_pro <- EFA(DOSPERT$cormat, n_factors = 10, N = DOSPERT$N,
                     type = "psych", rotation = "promax")
# EFAtools::EFA with type = "SPSS" with promax rotation
EFA_SPSS_pro <- EFA(DOSPERT$cormat, n_factors = 10, N = DOSPERT$N,
                    type = "SPSS", rotation = "promax")
# psych::fa with promax rotation
psych_pro <- psych::fa(DOSPERT$cormat, nfactors = 10, n.obs = DOSPERT$N,
                       fm = "pa", rotate = "Promax")

## Check replication of results

# Compare loadings from psych::fa and EFAtools::EFA with type = "psych"
COMPARE(EFA_psych_pro$rot_loadings, psych_pro$loadings)

# Compare loadings from SPSS and EFAtools::EFA with type = "SPSS"
COMPARE(EFA_SPSS_pro$rot_loadings, SPSS$DOSPERT$pro_load)

## Compare original results (just to see the difference)

# Compare loadings from psych::fa and SPSS
COMPARE(psych_pro$loadings, SPSS$DOSPERT$pro_load)


