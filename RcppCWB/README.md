
[![License: GPL
v3](http://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/RcppCWB)](https://cran.r-project.org/package=RcppCWB)
[![Travis-CI Build
Status](https://api.travis-ci.org/PolMine/RcppCWB.svg?branch=master)](https://travis-ci.org/PolMine/RcppCWB)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/PolMine/RcppCWB?branch=master&svg=true)](https://ci.appveyor.com/project/PolMine/RcppCWB)
[![codecov](https://codecov.io/gh/PolMine/RcppCWB/branch/master/graph/badge.svg)](https://codecov.io/gh/PolMine/RcppCWB/branch/dev)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3904235.svg)](https://doi.org/10.5281/zenodo.3904235)

# Rcpp bindings for the Corpus Workbench (CWB)

The package exposes functions of the Corpus Worbench (CWB) by way of
Rcpp wrappers. Furthermore, the packages includes Rcpp code for
performance critical operations. The main purpose of the package is to
serve as an interface to the CWB for the package
[polmineR](https://CRAN.R-project.org/package=RcppCWB).

There is a huge intellectual debt to the developers of the R-package
‘rcqp’, Bernard Desgraupes and Sylvain Loiseau. The main impetus for
developing RcppCWB is that using Rcpp decreases the pains to maintain
the package, to expand the CWB functionality exposed, and – most
importantly – to make it portable to Windows systems.

### Installation on Windows

Pre-compiled binaries of the package ‘RcppCWB’ can be obtained from
CRAN.

``` r
install.packages("RcppCWB")
```

If you want to get the development version, you need to compile RcppCWB
yourself. Having
[Rtools](https://cran.r-project.org/bin/windows/Rtools/) installed on
your system is necessary. Using the mechanism offered by the devtools
package, you can install RcppCWB from GitHub.

``` r
if (!"devtools" %in% installed.packages()[,"Package"]) install.packages("devtools")
devtools::install_github("PolMine/RcppCWB")
```

During the installation, cross-compiled versions of the corpus library
(CL) are downloaded from the GitHub repository
[PolMine/libcl](https://github.com/PolMine/libcl). It is not necessary
to install dependencies.

## Installation on Ubuntu

The package includes the source code of the Corpus Workbench (CWB),
slightly modified to make it compatible with R requirements. Compiling
the CWB requires the pcre and glib libraries to be present. Using the
Aptitude package manager (Ubuntu/Debian), running the following command
from the shell will fulfill these dependencies.

``` sh
sudo apt-get install libpcre3-dev libglib2.0-dev
```

Then, use the conventional R installation mechanism to install R
dependencies, and the release of RcppCWB at CRAN.

``` r
install.packages(pkgs = c("Rcpp", "knitr", "testthat"))
install.packages("RcppCWB")
```

To install the development version, using the mechanism offered by the
devtools package is recommended.

``` r
if (!"devtools" %in% installed.packages()[,"Package"]) install.packages("devtools")
devtools::install_github("PolMine/RcppCWB")
```

## Installation on MacOS

On macOS, the [pcre](http://www.pcre.org/) and
[Glib](https://developer.gnome.org/glib/) libraries need to be present.
We recommend to use ‘Homebrew’ as a package manager for macOS. To
install Homebrew, follow the instructions on the [Homebrew
Website](https://brew.sh/index_de.html). It may also be necessary to
also install [Xcode](https://developer.apple.com/xcode/) and
[XQuartz](https://www.xquartz.org).

The following commands then need to be executed from a terminal window.
They will install the C libraries the CWB relies on:

``` sh
brew -v install pkg-config
brew -v install glib --universal
brew -v install pcre --universal
brew -v install readline
```

Then open R and use the conventional R installation mechanism to install
dependencies, and the release of RcppCWB at CRAN.

``` r
install.packages(pkgs = c("Rcpp", "knitr", "testthat"))
install.packages("RcppCWB")
```

To install the development version, using the mechanism offered by the
devtools package is recommended.

``` r
if (!"devtools" %in% installed.packages()[,"Package"]) install.packages("devtools")
devtools::install_github("PolMine/RcppCWB")
```

## Usage

The package offers low-level access to CWB-indexed corpora. Using
RcppCWB may not intuitive at the outset: It is designed to serve as a an
efficient backend for packages offering higher-level functionality, such
as polmineR. the

RcppCWB includes a small sample corpus called (‘REUTERS’). After loading
the package, we need to determine whether we can use the registry
describing the corpus within the package, or whether we need to work
with a temporary registry.

``` r
library(RcppCWB)
if (!check_pkg_registry_files()){
  registry <- use_tmp_registry()
} else {
  registry <- get_pkg_registry()
} 
```

To start with, we get the number of tokens of the corpus.

``` r
cpos_total <- cl_attribute_size(
  corpus = "REUTERS", attribute = "word",
  attribute_type = "p", registry = registry
)
cpos_total
```

    ## [1] 4050

To decode the token stream of the corpus.

``` r
token_stream_str <- cl_cpos2str(
  corpus = "REUTERS", p_attribute = "word",
  cpos = seq.int(from = 0, to = cpos_total - 1),
  registry = registry
  )
```

To get the corpus positions of a token.

``` r
token_to_get <- "oil"
id_oil <- cl_str2id(corpus = "REUTERS", p_attribute = "word", str = token_to_get)
cpos_oil <- cl_id2cpos <- cl_id2cpos(corpus = "REUTERS", p_attribute = "word", id = id_oil)
```

Get the frequency of token.

``` r
oil_freq <- cl_id2freq(corpus = "REUTERS", p_attribute = "word", id = id_oil)
```

Using regular expressions.

``` r
ids <- cl_regex2id(corpus = "REUTERS", p_attribute = "word", regex = "M.*")
m_words <- cl_id2str(corpus = "REUTERS", p_attribute = "word", id = ids)
```

To use the CQP syntax, we need to initialize CQP first.

``` r
cqp_initialize(registry = registry)
```

    ## Warning in cqp_initialize(registry = registry): CQP has already been
    ## initialized. Re-initialization is not possible. Only resetting registry.

    ## [1] TRUE

``` r
cqp_query(corpus = "REUTERS", query = '"crude" "oil"')
```

    ## NULL

``` r
cpos <- cqp_dump_subcorpus(corpus = "REUTERS")
cpos
```

    ##       [,1] [,2]
    ##  [1,]   14   15
    ##  [2,]   56   57
    ##  [3,]  548  549
    ##  [4,]  584  585
    ##  [5,]  607  608
    ##  [6,] 2497 2498
    ##  [7,] 2842 2843
    ##  [8,] 2891 2892
    ##  [9,] 2928 2929
    ## [10,] 3644 3645
    ## [11,] 3709 3710
    ## [12,] 3998 3999

## License

The packge is licensed under the [GNU General Public
License 3](https://www.gnu.org/licenses/gpl-3.0.de.html). For the
copyrights for the ‘Corpus Workbench’ (CWB) and acknowledgement of
authorship, see the file COPYRIGHTS.

## Acknowledgements

There is a huge intellectual debt to the developers of the R-package
‘rcqp’, Bernard Desgraupes and Sylvain Loiseau. Developing RcppCWB
would have been unthinkable without their original work to wrap the CWB
into an R package.

The CWB is a classic and mature tool: The work of the CWB developers,
Oliver Christ, Bruno Maximilian Schulze, Arne Fitschen and Stefan Evert
is gratefully acknowledged.
