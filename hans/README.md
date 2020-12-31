# hans

<!-- badges: start -->
<!-- badges: end -->

The goal of hans is to make a fast dataframe compatible haversine function.

## Installation

You can install the released version of hans from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("hans")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(tibble)
library(dplyr)
library(magrittr)

set.seed(42)
lon1 <- runif(-160, -60, n = 1e6)
lat1 <- runif(40, 60, n = 1e6)
lon2 <- runif(-160, -60, n = 1e6)
lat2 <- runif(40, 60, n = 1e6)

df <- tibble::tibble(lat1, lon1, lat2, lon2) %>% 
  mutate(hav = haversine(lat1, lon1, lat2, lon2))
```

