
[![Project Status: Active ? The project has reached a stable, usable
state and is being actively
developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Build
Status](https://travis-ci.org/paulhibbing/AGread.svg?branch=master)](https://travis-ci.org/paulhibbing/AGread)
[![License](https://img.shields.io/badge/licence-MIT-blue.svg)](https://opensource.org/licenses/MIT)

-----

[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/AGread)](https://cran.r-project.org/package=AGread)

-----

<!-- README.md is generated from README.Rmd. Please edit that file -->

# AGread

AGread is for bringing ActiGraph sensor data into R. It is designed to
streamline and standardize the file-reading process, regardless of which
file format is provided (currently supporting `.gt3x` and `.csv`).

AGread can be used flexibly to develop new methods for handling
ActiGraph data, or to invoke existing methods, many of which exist in
other R packages
(e.g. [PhysicalActivity](https://cran.r-project.org/package=PhysicalActivity)
and [TwoRegression](https://cran.r-project.org/package=TwoRegression)),
or will soon.

New in AGread 1.0.0, `Rcpp` has been invoked to speed up the process of
reading `.gt3x` files. There is now documented equivalence between the
outcomes of `read_gt3x` and csv reading functions `read_AG__raw` and
`read_AG_IMU`.

## Installation

You can install the development version of AGread from github with:

``` r
# install.packages("devtools")
devtools::install_github("paulhibbing/AGread")
```

Windows users, make sure you have
[Rtools](https://cran.r-project.org/bin/windows/Rtools/) installed
before running the above.

Alternatively, AGread v1.0.0 is available on CRAN. Install it with:

``` r
install.packages("AGread")
```
