<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/cort)](https://CRAN.R-project.org/package=cort)
[![CRAN number of downloads](https://cranlogs.r-pkg.org/badges/grand-total/cort)](https://cranlogs.r-pkg.org/badges/grand-total/cort)
[![Travis build status](https://img.shields.io/travis/com/lrnv/cort/master?logo=travis&style=flat-square&label=Linux)](https://travis-ci.com/lrnv/cort)
[![Github Build status](https://img.shields.io/github/workflow/status/lrnv/cort/R%20CMD%20Check%20via%20%7Btic%7D?logo=github&label=Github%20build&style=flat-square)](https://github.com/lrnv/cort/actions)
[![AppVeyor build status](https://img.shields.io/appveyor/ci/lrnv/cort?label=Windows&logo=appveyor&style=flat-square)](https://ci.appveyor.com/project/lrnv/cort)
[![Codecov test coverage](https://codecov.io/gh/lrnv/cort/branch/master/graph/badge.svg)](https://codecov.io/gh/lrnv/cort?branch=master)
<!-- badges: end -->

The `cort` package provides S4 classes and methods to fit several copula models: 

* The classic empirical checkerboard copula and the empirical checkerboard copula with known margins, see Cuberos, Masiello and Maume-Deschamps (2019) are proposed. These two models allow to fit copulas in high dimension with a small number of observations, and they are always proper copulas. Some flexibility is added via a possibility to differentiate the checkerboard parameter by dimension. 

* The last model consist of the implementation of the Copula Recursive Tree algorithm, aka. CORT, including the localised dimension reduction, which fits a copula by recursive splitting of the copula domain, see Laverny, Maume-Deschamps, Masiello and Rullière (2020).

* We finally provide an efficient way of mixing copulas, allowing to bag the algorithm into a forest, and a generic way of measuring d-dimensional boxes with a given copula.

## Installation

`cort` is Now on [CRAN](https://CRAN.R-project.org)! You can install the stable version with:

``` r
install.packages("cort")
```

The upstream development version can also be installed with :

``` r
devtools::install_github("lrnv/cort")
```

Note that the installation from github will require the system to have a compiler: 

- Windows: Rtools
- macOS: Xcode CLI
- Linux: r-base-dev (debian)


The vignettes are quite expressive. They give a clear overview of what can be done with this package, how it is coded and why it is useful. Please read them for more details. 

## How to report bugs and get support

To report a bug, feel free to open an issue on the github repository. Support can also be provided through the same chanel if you need it.

## How to contribute

Every contribution is welcome, on the form of pull requests on the github repository. For large modifications, please open an issue for discussions firsts. Concerning the naming convention, the CamelCase functions usually designate classes and constructors of these classes, and all other methods are in snake_case.


## References

Cuberos A, Masiello E, Maume-Deschamps V (2019). “Copulas Checker-Type Approximations: Application to Quantiles Estimation of Sums of Dependent Random Variables.” *Communications in Statistics - Theory and Methods, 1--19. ISSN 0361-0926, 1532-415X.*

Laverny O, Maume-Deschamps V, Masiello E, Rullière D (2020). “Dependence Structure Estimation Using Copula Recursive Trees.” *arXiv preprint arXiv:2005.02912*
