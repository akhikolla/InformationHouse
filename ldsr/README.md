<!-- badges: start -->
[![Travis build status](https://travis-ci.org/ntthung/ldsr.svg?branch=master)](https://travis-ci.org/ntthung/ldsr)
[![GPL license](https://img.shields.io/badge/License-GPL-blue.svg)](http://perso.crans.org/besson/LICENSE.html)
[![GitHub version](https://badge.fury.io/gh/ntthung%2Fldsr.svg)](https://badge.fury.io/gh/ntthung%2Fldsr)
<!-- badges: end -->


# ldsr

Streamflow reconstruction using linear dynamical system. The method was described in Nguyen and Galelli (2018) and used in Nguyen et al (2019). Please cite Nguyen and Galelli (2018) when using the package.

To install this package, type the following commands into your RStudio console:

```
install.packages('devtools')
devtools::install_github('ntthung/ldsr', build_vignettes = TRUE)
```

Note that if you don't already have Rtools (a set of tools to build packages from sources), RStudio will prompt you to install it. Just follow the instructions. After Rtools installation is complete, type the above `install_github` command again.

To view the package's vignette, type

`browseVignettes('ldsr')`

# References

Nguyen, H. T. T., & Galelli, S. (2018). A linear dynamical systems approach to streamflow reconstruction reveals history of regime shifts in northern Thailand. Water Resources Research, 54, 2057– 2077. https://doi.org/10.1002/2017WR022114 

Nguyen, H. T. T., Turner, S. W., Buckley, B. M., & Galelli, S. (2019). Coherent streamflow variability in Monsoon Asia over the past eight centuries---links to oceanic drivers. https://doi.org/10.31223/osf.io/5tg68
