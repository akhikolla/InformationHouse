
<!-- README.md is generated from README.Rmd. Please edit that file -->

# triangulr <a href="https://irkaal.github.io/triangulr/"><img src="man/figures/logo.png" align="right" height="139" /></a>

[![Build](https://github.com/irkaal/triangulr/workflows/R-CMD-check/badge.svg)](https://github.com/irkaal/triangulr/actions)
[![Codecov](https://codecov.io/gh/irkaal/triangulr/branch/master/graph/badge.svg)](https://codecov.io/gh/irkaal/triangulr?branch=master)
[![CRAN](http://www.r-pkg.org/badges/version/triangulr)](https://cran.r-project.org/package=triangulr)
[![Downloads](http://cranlogs.r-pkg.org/badges/triangulr)](https://cran.r-project.org/package=triangulr)
[![Downloads
Overall](http://cranlogs.r-pkg.org/badges/grand-total/triangulr)](https://cran.r-project.org/package=triangulr)
[![Lifecycle](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://www.tidyverse.org/lifecycle/#stable)

## Introduction

The `triangulr` package provides high-performance triangular
distribution functions which includes density function, distribution
function, quantile function, random variate generator, moment generating
function, characteristic function, and expected shortfall function for
the triangular distribution.

## Installation

You can install the released version of `triangulr` from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("triangulr")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("irkaal/triangulr")
```

## Example

These are basic examples of using the included functions:

``` r
library(triangulr)
```

Using the density function,
[`dtri()`](https://irkaal.github.io/triangulr/reference/Triangular.html).

``` r
x <- c(0.1, 0.5, 0.9)

dtri(x,
     min = 0,
     max = 1,
     mode = 0.5)
#> [1] 0.4 2.0 0.4

dtri(x,
     min = c(0, 0, 0),
     max = 1,
     mode = 0.5)
#> [1] 0.4 2.0 0.4
```

Using the distribution function,
[`ptri()`](https://irkaal.github.io/triangulr/reference/Triangular.html).

``` r
q <- c(0.1, 0.5, 0.9)

1 - ptri(q, lower_tail = FALSE)
#> [1] 0.02 0.50 0.98

ptri(q, lower_tail = TRUE)
#> [1] 0.02 0.50 0.98

ptri(q, log_p = TRUE)
#> [1] -3.91202301 -0.69314718 -0.02020271

log(ptri(q, log_p = FALSE))
#> [1] -3.91202301 -0.69314718 -0.02020271
```

Using the quantile function,
[`qtri()`](https://irkaal.github.io/triangulr/reference/Triangular.html).

``` r
p <- c(0.1, 0.5, 0.9)

qtri(1 - p, lower_tail = FALSE)
#> [1] 0.2236068 0.5000000 0.7763932

qtri(p, lower_tail = TRUE)
#> [1] 0.2236068 0.5000000 0.7763932

qtri(log(p), log_p = TRUE)
#> [1] 0.2236068 0.5000000 0.7763932

qtri(p, log_p = FALSE)
#> [1] 0.2236068 0.5000000 0.7763932
```

Using the random variate generator,
[`rtri()`](https://irkaal.github.io/triangulr/reference/Triangular.html).

``` r
n <- 3

set.seed(1)
rtri(n,
     min = 0,
     max = 1,
     mode = 0.5)
#> [1] 0.3643547 0.4313490 0.5378601

set.seed(1)
rtri(n,
     min = c(0, 0, 0),
     max = 1,
     mode = 0.5)
#> [1] 0.3643547 0.4313490 0.5378601

dqrng::dqset.seed(1)
rtri(
  n,
  min = c(0, 0, 0),
  max = 1,
  mode = 0.5,
  dqrng = TRUE
)
#> [1] 0.3951856 0.8516496 0.4494472
```

Using the moment generating function,
[`mgtri()`](https://irkaal.github.io/triangulr/reference/Triangular.html).

``` r
t <- c(1, 2, 3)

mgtri(t,
      min = 0,
      max = 1,
      mode = 0.5)
#> [1] 1.683357 2.952492 5.387626

mgtri(t,
      min = c(0, 0, 0),
      max = 1,
      mode = 0.5)
#> [1] 1.683357 2.952492 5.387626
```

Using the characteristic function,
[`ctri()`](https://irkaal.github.io/triangulr/reference/Triangular.html).

``` r
t <- c(1, 2, 3)

ctri(t,
     min = 0,
     max = 1,
     mode = 0.5)
#> [1] 0.8594513+0.4695204i 0.4967514+0.7736445i 0.0584297+0.8239422i

ctri(t,
     min = c(0, 0, 0),
     max = 1,
     mode = 0.5)
#> [1] 0.8594513+0.4695204i 0.4967514+0.7736445i 0.0584297+0.8239422i
```

Using the expected shortfall function,
[`estri()`](https://irkaal.github.io/triangulr/reference/Triangular.html).

``` r
p <- c(0.1, 0.5, 0.9)

estri(p,
      min = 0,
      max = 1,
      mode = 0.5)
#> [1] 0.1490712 0.3333333 0.4610079

estri(p,
      min = c(0, 0, 0),
      max = 1,
      mode = 0.5)
#> [1] 0.1490712 0.3333333 0.4610079
```
