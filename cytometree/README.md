
<!-- README.Rmd generates from README.Rmd. Please generate that file once done editing -->

# cytometree

[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/cytometree)](https://cran.r-project.org/package=cytometree)
[![Travis-CI Build
Status](https://travis-ci.org/sistm/cytometree.svg?branch=master)](https://travis-ci.org/sistm/cytometree)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/sistm/cytometree?branch=master&svg=true)](https://ci.appveyor.com/project/borishejblum/cytometree)
[![codecov.io](https://codecov.io/github/sistm/Cytometree/coverage.svg?branch=master)](https://codecov.io/github/sistm/Cytometree?branch=master)
[![Downloads](https://cranlogs.r-pkg.org/badges/cytometree?color=blue)](https://www.r-pkg.org/pkg/cytometree)

## Overview

`cytometree` is a package which performs **automatic gating and
annotation of flow-cytometry data**. On top of the [CRAN help
files](https://cran.r-project.org/package=cytometree/cytometree.pdf), we
also provide a
[vignette](https://cran.r-project.org/package=cytometree/vignettes/autogating_cytometree.html)
illustrating the functionalities of `cytometree`.

The `cytometree` algorithm rely on the construction of a **binary
tree**, the nodes of which represents **cellular (sub)populations**. At
each node, observed cellular markers are modeled by both a family of
normal and a family of normal mixture distributions and splitting of
cells into further subpopulations is decided according to a normalized
difference of AIC between the two families.

Given the **unsupervised** nature of such a binary tree, some of the
available markers may not be used to find the different cell populations
present in a given sample. So in order to recover a complete annotation,
we propose a **post processing annotation** procedure which allows the
user to distinguish two or three expression levels per marker.

The following article explains in more details how `cytometree` works:

> Commenges D, Alkhassim C, Gottardo R, Hejblum BP, Thiébaut R (2018).
> cytometree: a binary tree algorithm for automatic gating in cytometry
> analysis. *Cytometry Part A* **93**(11):1132-1140.
> [\<doi: 10.1002/cyto.a.23601\>](https://doi.org/10.1002/cyto.a.23601)

## Installation

The easiest way to get `cytometree` is to install it from
[CRAN](https://cran.r-project.org/package=cytometree):

``` r
install.packages("cytometree")
```

Or to get the development version from
[GitHub](https://github.com/sistm/cytometree):

``` r
#install.packages("devtools")
devtools::install_github("sistm/cytometree")
```

– Chariff Alkhassim & Boris Hejblum
