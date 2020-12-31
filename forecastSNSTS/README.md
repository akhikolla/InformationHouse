[![Build Status](https://travis-ci.org/tobiaskley/forecastSNSTS.svg?branch=develop)](https://travis-ci.org/tobiaskley/forecastSNSTS)
[![Coverage Status](https://img.shields.io/codecov/c/github/tobiaskley/forecastSNSTS/develop.svg)](https://codecov.io/github/tobiaskley/forecastSNSTS?branch=develop)
[![CRAN RStudio mirror downloads](http://cranlogs.r-pkg.org/badges/forecastSNSTS)](https://cran.r-project.org/package=forecastSNSTS)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/forecastSNSTS)](https://cran.r-project.org/package=forecastSNSTS)


forecastSNSTS: Forecasting of Stationary and Non-Stationary Time Series
=======================================================================

The `forecastSNSTS` package provides methods to compute linear h-step prediction coefficients based on localised and iterated Yule-Walker estimates and empirical mean square prediction errors from the resulting predictors.

It is intended to support the paper
[Predictive, finite-sample model choice for time series under stationarity and non-stationarity](https://arxiv.org/abs/1611.04460), which we refer to as Kley et al. (2019).

You can track (and contribute to) the development of `forecastSNSTS` at https://github.com/tobiaskley/forecastSNSTS. If you encounter unexpected behaviour while using `forecastSNSTS`, please write an [email](mailto:tobias.kley@bristol.ac.uk) or file an [issue](http://github.com/tobiaskley/forecastSNSTS/issues).

## Getting started with ``forecastSNSTS``

First, if you have not done so already, install R from http://www.r-project.org (click on download R, select a location close to you, and download R for your platform). Once you have the latest version of R installed and started execute the following commands on the R shell:

 ```
 install.packages("forecastSNSTS")
 devtools::install_github("tobiaskley/forecastSNSTS", ref="develop")
 ```

This will first install the R package ``devtools`` and then use it to install the latest (development) version of ``forecastSNSTS`` from the GitHub repository. In case you do not have LaTeX installed on your computer you may want to use

Now that you have R and ``forecastSNSTS`` installed you can access all the functions available. To load the package and access the help files:

```
library(forecastSNSTS)
help("forecastSNSTS")
```

A demo is available. It can be started by

```
demo("tvARMA11")
```

At the bottom of the online help page to the package you will find an index to all the help files available.

## Replicating the examples of the paper with ``forecastSNSTSexamples``

Note that there is a separate R package, called [forecastSNSTSexamples](http://github.com/tobiaskley/forecastSNSTSexamples) and available only on GitHub, that can be used to replicate the empirical examples from Section 5 of Kley et al. (2019). 
