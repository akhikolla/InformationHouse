[![Build Status](https://travis-ci.org/katzfuss-group/GPvecchia.svg?branch=master)](https://travis-ci.org/katzfuss-group/GPvecchia)
[![CRAN status](https://www.r-pkg.org/badges/version/tidyverse)](https://cran.r-project.org/package=tidyverse)

# GPvecchia
Fast Gaussian-process inference using general Vecchia approximations

For examples of how to use the package, please see the vignettes folder. Please note that GPvecchia is under active development and not stable at this time.

## Reporting problems
If you have an issue with GPvecchia, including unsatisfactory runtime, please open a Github ticket.

## References
[Katzfuss, M., & Guinness, J. (2017). A general framework for Vecchia approximations of Gaussian processes. *arXiv:1708.06302*.](https://arxiv.org/abs/1708.06302)

[Katzfuss, M., Guinness, J., Gong, W., & Zilber, D. (2018). Vecchia approximations of Gaussian-process predictions. *arXiv:1805.03309*.](https://arxiv.org/abs/1805.03309)

[Zilber, D., & Katzfuss, M. (2019). Vecchia-Laplace approximations of generalized Gaussian processes for big non-Gaussian spatial data. *arXiv:1906.07828*.](https://arxiv.org/abs/1906.07828)

## Installation
<!---To ensure that the algorithms run efficiently on your computer, we recommend installing the package by downloading the repo to a local folder, setting your R working directory to that folder, and then running the following code in R:
```{r}
library(GpGp); library(Matrix); library(RcppParallel)
library(parallel); library(sparseinv); library(fields)
for (nm in list.files('GPvecchia/R',pattern = "\\.[RrSsQq]$")) {
  cat(nm,":"); source(file.path('GPvecchia/R',nm)); cat("\n")
}
Rcpp::sourceCpp('GPvecchia/src/U_NZentries.cpp')
Rcpp::sourceCpp('GPvecchia/src/MaxMin.cpp')
```
--->
 
This package can be installed directly from Github by running
```{r}
devtools::install_github("katzfuss-group/GPvecchia")
```
Alternatively, one can download the repository and then build the package manually:
```{bash}
R CMD build GPvecchia
R CMD INSTALL GPvecchia_0.1_0.tar.gz
```

<!---
.tar.gz file from the main directory here and then run:
```{r}
install.packages("GPvecchia_0.1.0.tar.gz", repos = NULL, type = "source")
```
-->

Note that [Rtools](https://cran.r-project.org/bin/windows/Rtools/) is required for compiling C/C++ with OpenMP on Windows systems. When installing Rtools, the system PATH needs to be set so that the C++ compiler included in Rtools can be found by R. Once Rtools is installed, `system('g++ -v')` can be used to check if the compiler is accessible from R.

