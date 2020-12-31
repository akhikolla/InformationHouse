# serrsBayes 0.4-1

## New Features

* Added example of `fitSpectraSMC` with informative priors to the methanol vignette.

## Bug Fixes

* Removed the R package ``hyperSpec`` from Suggests because it is no longer available on CRAN.
* Fixed dependency on ``hyperSpec`` in Introduction vignette due to new CRAN check. This also affects the ``TAMRA`` dataset.
* Fixed a rare bug in reweighting that could potentially result in an infinite loop.
* Peak locations are now constrained to lie within the range of observed wavelengths.
* Locations are sorted to preserve identifiability.
* Fixed error in equations for FWHM in `mixedVoigt` & `getVoigtParam` as well as in both vignettes.

# serrsBayes 0.4-0

## New Features

* OpenMP implementation of mhUpdateVoigt provides around 3x improvement in elapsed runtime
* New vignette illustrates when to use the 3 functions `fitSpectraMCMC`, `fitSpectraSMC`, and `fitVoigtPeaksSMC`
* New datasets ``methanol`` and ``TAMRA`` for use in the vignettes
* Debug output prints the means of the amplitudes of the peaks at each SMC iteration

## Bug Fixes

* Fixed bug in ``fitVoigtPeaksSMC`` where it was using the upper-triangular instead of lower-triangular Cholesky factorisation to generate random-walk Metropolis proposals
* Introduction vignette no longer attempts to download a .zip file over HTTPS
* Deleted unused Rcpp function randomWalkVoigt

# serrsBayes 0.3-13

* changed maintainer email address due to new academic affiliation
* hex sticker in README

# serrsBayes 0.3-12

## New Features

* Vignette with example for tetramethylrhodamine (TAMRA)
* pkgdown website https://mooresm.github.io/serrsBayes/

## Bug Fixes

* Fixed compile errors on Solaris

# serrsBayes 0.3-10

* First version released on CRAN
