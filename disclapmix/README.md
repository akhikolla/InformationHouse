# disclapmix

An R package (<https://www.r-project.org/>) to perform inference in a mixture of discrete Laplace distributions using the EM algorithm.
Intended for forensic Y chromosomal STR (Y-STR) haplotype analyses. 

See documentation included in package (vignettes and manual) at <https://mikldk.github.io/disclapmix/>.


## Installation

You first need `R` (<https://www.r-project.org/>). 
Then you can install `disclapmix` from CRAN using

```r
install.packages("disclapmix")
```

You can also install the development version from GitHub by using the `remotes` package (<https://CRAN.R-project.org/package=remotes>):

```r
# install.packages("remotes")
remotes::install_github("mikldk/disclapmix")
```

## Getting started

Refer to the included vignettes. You can get an overview of the included vignettes by the following `R` command:

```r
vignette(package = "disclapmix")
```

To read a vignette, type:

```r
vignette("introduction", package = "disclapmix")
```

### Running tests

Note that to also install the tests, you need to install the package as follows:

``` r
# install.packages("remotes")
remotes::install_github("mikldk/disclapmix", INSTALL_opts="--install-tests")
```

You can now run the tests:

``` r
library(disclapmix)
library(testthat)
test_package('disclapmix')
```

## Contribute, issues, and support

Please use the issue tracker at <https://github.com/mikldk/disclapmix/issues> 
if you want to notify us of an issue or need support.
If you want to contribute, please either create an issue or make a pull request.


## References

Andersen MM, PS Eriksen, N Morling (2013). 
*The discrete Laplace exponential family and estimation of Y-STR haplotype frequencies*.
Journal of Theoretical Biology 329. 
<https://doi.org/10.1016/j.jtbi.2013.03.009>

Andersen MM, PS Eriksen, N Morling (2013). 
*A gentle introduction to the discrete Laplace method for estimating Y-STR haplotype frequencies*.
arXiv:1304.2129. 
<https://arxiv.org/abs/1304.2129>

## Disclaimer

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

## License

License: GPL-2.

## Badges

The Journal of Open Source Software:

[![DOI](http://joss.theoj.org/papers/10.21105/joss.00748/status.svg)](https://doi.org/10.21105/joss.00748)

Zenodo: 

[![DOI](https://zenodo.org/badge/130306482.svg)](https://zenodo.org/badge/latestdoi/130306482)

Travis CI:

[![Travis-CI Build Status](https://travis-ci.org/mikldk/disclapmix.svg?branch=master)](https://travis-ci.org/mikldk/disclapmix)

