# easyVerification

This package provides functions to simplify application of forecast verification metrics to large datasets of ensemble forecasts. The design goals of `easyVerification` are:

* __Flexibility:__ a variety of data structures are supported
* __Ease of use:__ Absolute forecasts and observations are converted to category and probability forecasts based on the threshold or probability (e.g. terciles) provided, ouputs are reformatted to fit the input
* __Convenience and flexibility over speed:__ R's built-in vectorisation is used where possible but more importantly, new metrics should be easy to implement

The forecast metrics are imported from the `SpecsVerification` package. Additional verification metrics not available through `SpecsVerification` are implemented directly. At the time of publication, the package offers functionality to compute the following deterministic and probabilitistic scores and skill scores:

1. Mean error (`EnsMe`), mean absolute error(`EnsMae`), mean squared error (`EnsMse`), and root mean squared error (`EnsRmse`) of the ensemble mean and their skill scores (e.g. `EnsRmsess`)
2. Correlation with the ensemble mean (`EnsCorr`)
3. Spread to error ratio (`EnsSprErr` and `FairSprErr`)
4. Area under the ROC curve (`EnsRoca`) and its skill score (`EnsRocss`)
5. Fair (`FairRps`) and standard (`EnsRps`) rank probability scores and skill scores (e.g. `FairRpss`)
6. Fair (`FairCrps`) and standard (`EnsCrps`) continuous ranked probability scores and skill scores (e.g. `FairCrpss`)
7. Dressed scores (`DressIgn`, `DressCrps`) and their skill scores (`DressIgnSs`, `DressCrpss`) with default ensemble dressing method ("silverman")
8. The generalized discrimination score for ensembles (`Ens2AFC`)

Additional forecast verification metrics can be added by the user following the examples above.

## Installation 
You can get the latest version from CRAN

```r
install.packages("easyVerification")
```

You can get the latest development version using

```r
devtools::install_github("MeteoSwiss/easyVerification")
```

## Getting started

You can find out more about the package and its functionality in the vignette.


```r
vignette('easyVerification')
```

The following example illustrates how to compute the continous ranked probability skill score of an ensemble forecast:


```r
suppressPackageStartupMessages(library(easyVerification))

## check out what is included in easyVerification
ls(pos="package:easyVerification")
#>  [1] "convert2prob" "count2prob"   "Ens2AFC"      "EnsCorr"     
#>  [5] "EnsError"     "EnsErrorss"   "EnsMae"       "EnsMaess"    
#>  [9] "EnsMe"        "EnsMess"      "EnsMse"       "EnsMsess"    
#> [13] "EnsRmse"      "EnsRmsess"    "EnsRoca"      "EnsRocss"    
#> [17] "EnsSprErr"    "FairSprErr"   "toyarray"     "toymodel"    
#> [21] "veriApply"

## set up the forecast and observation data structures
## assumption: we have 13 x 5 spatial instances, 15 forecast 
## times and 51 ensemble members
tm <- toyarray(c(13,5), N=15, nens=51)
fo.crpss <- veriApply("EnsCrpss", fcst=tm$fcst, obs=tm$obs)

## if the data are organized differently such that forecast
## instance and ensemble members are NOT the last two array
## dimensions, this has to be indicated

## alternative setup:
## forecast instance, ensemble members, all forecast locations
## collated in one dimension
fcst2 <- array(aperm(tm$fcst, c(3,4,1,2)), c(15, 51, 13*5))
obs2 <- array(aperm(tm$obs, c(3,1,2)), c(15, 13*5))
fo2.crpss <- veriApply("EnsCrpss", fcst=fcst2, obs=obs2, 
                       ensdim=2, tdim=1)

## The forecast evaluation metrics are the same, but the 
## data structure is different in the two cases
dim(fo.crpss$crpss)
#> [1] 13  5
dim(fo2.crpss$crpss)
#> [1] 65
range(fo.crpss$crpss - c(fo2.crpss$crpss))
#> [1] 0 0
```

## Parallel processing

As of `easyVerification 0.1.7.0`, parallel processing is supported under *NIX systems. The following minimal example illustrates how to use the parallel processing capabilities of `easyVerification`.

```r
## generate a toy-model forecast observation set of 
## 10 x 10 forecast locations (e.g. lon x lat)
tm <- toyarray(c(10,10))

## run and time the ROC skill score for tercile forecasts without parallelization
system.time({
  tm.rocss <- veriApply("EnsRocss", tm$fcst, tm$obs, prob=1:2/3)  
})
#>    user  system elapsed 
#>   1.984   0.008   2.000

## run the ROC skill score with parallelization
system.time({
  tm.rocss.par <- veriApply("EnsRocss", tm$fcst, tm$obs, prob=1:2/3, parallel=TRUE)
})
#> Loading required namespace: parallel
#> [1] "Number of CPUs 3"
#>    user  system elapsed 
#>   0.088   0.040   0.824
```

To get additional help and examples please see the vignette `vignette('easyVerification')` or the help pages of the functions in `easyVerification` (e.g. `help(veriApply)`).
