
<!-- README.md is generated from README.Rmd. Please edit that file -->

# waspr

<!-- [![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/waspr)](https://cran.r-project.org/package=waspr) -->

[![Build
Status](https://travis-ci.org/joliencremers/waspr.svg?branch=master)](https://travis-ci.org/joliencremers/waspr)

The goal of waspr is to compute Wasserstein barycenters of subset
posteriors.

## Installation

The R-package waspr can be installed from CRAN as follows:

``` r
install.packages("waspr")
```

You can install a beta-version of waspr from github with:

``` r
install.packages("devtools")
devtools::install_github("joliencremers/waspr")
```

## Example

This is a basic example which shows you how to compute the Wasserstein
barycenter from a set of MCMC outputs for several data subsets. A more
extensive explanation of the usage of the package can be found in the
Tutorial vignette.

``` r
library(waspr)
#> 
#> Attaching package: 'waspr'
#> The following object is masked from 'package:base':
#> 
#>     summary

wasp(pois_logistic,
     par.names = c("beta_s", "alpha_l", "beta_l",
                   "baseline_sigma", "baseline_mu",
                   "correlation", "sigma_s", "sigma_l"))
#> 
#> 
#> WASP 
#> 
#> Call: 
#> wasp(mcmc = pois_logistic, par.names = c("beta_s", "alpha_l", 
#>     "beta_l", "baseline_sigma", "baseline_mu", "correlation", 
#>     "sigma_s", "sigma_l"))
#> 
#> Swapping algorithm: 
#> iter = 10
#> acc = 0.001
#> 
#> MCMC: 
#> subsets = 8
#> parameters = 8
#> samples = 450
#> 
#> Posterior summary of the Wasserstein Barycenter: 
#>                      mean       mode         sd      LB HPD     UB HPD
#> beta_s          0.5527601  0.5518034 0.10988949  0.36598187  0.7896041
#> alpha_l         2.6811079  2.6959176 0.19199304  2.30380675  3.0295802
#> beta_l          0.7508520  0.7339988 0.21631011  0.37281283  1.1740767
#> baseline_sigma  0.3563222  0.3811609 0.06859910  0.21910807  0.4870079
#> baseline_mu    -0.8008872 -0.7516167 0.10867533 -1.01168299 -0.5944583
#> correlation     0.1732170  0.1392670 0.07437737  0.02824474  0.3059979
#> sigma_s         1.7225455  1.7535499 0.17920847  1.40126462  2.0610585
#> sigma_l         1.2190297  1.2612822 0.07558163  1.06768047  1.3569757
```
