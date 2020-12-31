## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.align = "center"
)

if (!requireNamespace("microbenchmark", quietly = TRUE)) {
      stop("Package \"microbenchmark\" needed for this vignette to work. Please install it.",
      call. = FALSE)
}


## -----------------------------------------------------------------------------
library(EFAtools)

## -----------------------------------------------------------------------------
# only use a subset to make analyses faster
DOSPERT_sub <- DOSPERT_raw[1:500,]

## -----------------------------------------------------------------------------
# Bartlett's test of sphericity
BARTLETT(DOSPERT_sub)

# KMO criterion
KMO(DOSPERT_sub)

## -----------------------------------------------------------------------------
# determine the number of factors to retain using parallel analysis
PARALLEL(DOSPERT_sub, eigen_type = "SMC")

## -----------------------------------------------------------------------------
# determine the number of factors to retain using parallel analysis
print(PARALLEL(DOSPERT_sub, eigen_type = "SMC"), plot = FALSE)

## -----------------------------------------------------------------------------
# determine the number of factors to retain using parallel analysis
print(EKC(DOSPERT_sub), plot = FALSE)

## -----------------------------------------------------------------------------
N_FACTORS(DOSPERT_sub, criteria = c("PARALLEL", "EKC", "SMT"),
          eigen_type_other = c("SMC", "PCA"))

## -----------------------------------------------------------------------------
N_FACTORS(DOSPERT_sub, method = "ULS")

## -----------------------------------------------------------------------------
N_FACTORS(test_models$baseline$cormat, N = 500,
          method = "ULS", eigen_type_other = c("SMC", "PCA"))

## -----------------------------------------------------------------------------
EFA(DOSPERT_sub, n_factors = 6)

## -----------------------------------------------------------------------------
EFA(DOSPERT_sub, n_factors = 6, rotation = "promax")

## -----------------------------------------------------------------------------
EFA(DOSPERT_sub, n_factors = 6, rotation = "promax", type = "psych")

## -----------------------------------------------------------------------------
EFA(DOSPERT_sub, n_factors = 6, rotation = "promax", type = "SPSS")

## -----------------------------------------------------------------------------
COMPARE(
  EFA(DOSPERT_sub, n_factors = 6, rotation = "promax", type = "psych")$rot_loadings,
  EFA(DOSPERT_sub, n_factors = 6, rotation = "promax", type = "SPSS")$rot_loadings
)


## -----------------------------------------------------------------------------
EFA(DOSPERT_sub, n_factors = 6, rotation = "oblimin", method = "ULS")

## -----------------------------------------------------------------------------
COMPARE(
  EFA(DOSPERT_sub, n_factors = 6, rotation = "promax")$rot_loadings,
  EFA(DOSPERT_sub, n_factors = 6, rotation = "oblimin", method = "ULS")$rot_loadings,
  x_labels = c("PAF and promax", "ULS and oblimin")
)

## -----------------------------------------------------------------------------
EFA_mod <- EFA(DOSPERT_sub, n_factors = 6, rotation = "promax")
fac_scores <- FACTOR_SCORES(DOSPERT_sub, f = EFA_mod)

## ----message=FALSE------------------------------------------------------------
microbenchmark::microbenchmark(
  PARALLEL(DOSPERT_sub, eigen_type = "SMC", n_datasets = 25),
  psych::fa.parallel(DOSPERT_sub, SMC = TRUE, plot = FALSE, n.iter = 25)
)

## ----message=FALSE------------------------------------------------------------
microbenchmark::microbenchmark(
  EFA(DOSPERT_raw, 6),
  psych::fa(DOSPERT_raw, 6, rotate = "none", fm = "pa")
)

## -----------------------------------------------------------------------------
# Average solution across many different EFAs with oblique rotations
EFA_AV <- EFA_AVERAGE(test_models$baseline$cormat, n_factors = 3, N = 500,
                      method = c("PAF", "ML", "ULS"), rotation = "oblique",
                      show_progress = FALSE)

# look at solution
EFA_AV

## -----------------------------------------------------------------------------
efa_dospert <- EFA(DOSPERT_sub, n_factors = 6, rotation = "promax")
efa_dospert

## -----------------------------------------------------------------------------
sl_dospert <- SL(efa_dospert)
sl_dospert

## -----------------------------------------------------------------------------
OMEGA(sl_dospert, type = "psych")

## -----------------------------------------------------------------------------
OMEGA(sl_dospert, factor_corres = matrix(c(rep(0, 18), rep(1, 6), rep(0, 30), 
                                         rep(1, 6), rep(0, 6), 1, 0, 1, 0, 1,
                                         rep(0, 19), rep(1, 6), rep(0, 31), 1, 0,
                                         1, 0, 1, rep(0, 30), rep(1, 6), 
                                         rep(0, 12)), ncol = 6, byrow = FALSE))

