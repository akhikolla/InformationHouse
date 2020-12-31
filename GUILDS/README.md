# GUILDS
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/GUILDS)](https://cran.r-project.org/package=GUILDS)
[![Build Status](https://travis-ci.org/thijsjanzen/GUILDS.svg?branch=master)](https://travis-ci.org/thijsjanzen/GUILDS)
[![codecov](https://codecov.io/gh/thijsjanzen/GUILDS/branch/master/graph/badge.svg)](https://codecov.io/gh/thijsjanzen/GUILDS)
[![](http://cranlogs.r-pkg.org/badges/grand-total/GUILDS)](https://cran.r-project.org/package=GUILDS)
[![](http://cranlogs.r-pkg.org/badges/GUILDS)](https://cran.r-project.org/package=GUILDS)

The GUILDS package combines a range of sampling formulas for the unified neutral model of biogeography and biodiversity. Alongside the sampling formulas, it includes methods to perform maximum likelihood optimization of the sampling formulas, methods to generate data given the neutral model, and methods to estimate the expected species abundance distribution. Sampling formulas included in the GUILDS package are the Etienne Sampling Formula (Etienne 2005), the guild sampling formula, where guilds are assumed to differ in dispersal ability (Janzen et al. 2015), and the guilds sampling formula conditioned on guild size (Janzen et al. 2015).

Furthermore it contains functions to generate data given the guilds model, with or without conditioning on guild size. C++ Code to obtain Sterling numbers of the first kind was adopted from the Tetame program by Jabot et al. (2008). 

### Updates 
- Version 1.3
  - GUILDS is now on GitHub: https://github.com/thijsjanzen/GUILDS
  - Wrote code tests to check code integrity, code coverage is >95
  - Modified maximum likelihood functions to take into account theta_x = theta_y = theta / 2
  - added a function to plot preston-style plots
- Version 1.2.1
  - Updated the User manual
- Version 1.2
  - fixed memory leak issues by adding extra vector access checks
  - fixed memory leak issues by introducing vectors in KDA code
  - renamed logLik to avoid shadowing of the function logLik in the package stats
- Version 1.1
  - removed malloc header from KDA code

### References
Janzen, T., Haegeman B., Etienne, R.S. (2015) A sampling formula for communities with multiple dispersal syndromes. Journal of Theoretical Biology 374: 94-106

Etienne, R.S. (2005). A new sampling formula for neutral biodiversity. Ecology Letters, 8(3), 253-260.

Jabot, F., Etienne, R.S., & Chave, J. (2008). Reconciling neutral community models and environmental filtering: theory and an empirical test. Oikos 117: 1308-1320


