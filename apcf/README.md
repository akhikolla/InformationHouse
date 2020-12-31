## apcf: Adapted Pair Correlation Function

[![Travis-CI Build Status](https://travis-ci.org/rnuske/apcf.svg?branch=master)](https://travis-ci.org/rnuske/apcf) 
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/rnuske/apcf?branch=master&svg=true)](https://ci.appveyor.com/project/rnuske/apcf) 
[![Drone.io Status](https://cloud.drone.io/api/badges/rnuske/apcf/status.svg)](https://cloud.drone.io/rnuske/apcf) 
[![Package-License](https://img.shields.io/badge/license-GPL--3-brightgreen.svg?style=flat)](https://www.gnu.org/licenses/gpl-3.0.html) 
[![CRAN](https://www.r-pkg.org/badges/version/apcf)](https://cran.r-project.org/package=apcf) 
[![Dependencies](https://tinyverse.netlify.com/badge/apcf)](https://cran.r-project.org/package=apcf) 
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2535612.svg)](https://doi.org/10.5281/zenodo.2535612) 


The Adapted Pair Correlation Function transfers the concept of the Pair Correlation Function from point patterns to patterns of patches of finite size and irregular shape (eg. lakes within a country). The main tasks are (i) the construction of nullmodels by rondomizing the patches of the original pattern within the study area, (ii) the edge correction by determining the proportion of a buffer within the study area, and (iii) the calculation of the shortest distances between the patches.

This is a reimplementation of the Adapted Pair Correlation Function (Nuske et al. 2009) in C++ using the libraries GEOS and GDAL directly instead of through PostGIS.


### Requirements
For Unix-alikes GDAL (>= 2.0.1) and GEOS (>= 3.4.0) are required.

On Ubuntu bionic (18.04) and beyond one can install the dependencies simply with `sudo apt install libgdal-dev libgeos-dev`. 
In earlier Ubuntu version either add [ubuntugis-unstable](http://ppa.launchpad.net/ubuntugis/ubuntugis-unstable/ubuntu/) to the `sources.list` and use above command or compile dependencies from source.


### Installation
The stable version can be installed from CRAN
```r
install.packages("apcf")
```

and the development is available from Github using the package remotes (formerly devtools)
```r
if(!require("remotes")) install.packages("remotes")
remotes::install_github("rnuske/apcf")
```


### Usage
```r
# calculate distances between patches of original pattern and 3 nullmodels
# number of nullmodels should by at least 199 and better yet 999
ds <- pat2dists(area=system.file("shapes/sim_area.shp", package="apcf"),
                pattern=system.file("shapes/sim_pat_reg.shp", package="apcf"),
                max_dist=25, n_sim=3)

# derive PCF and envelope from distances
pcf <- dists2pcf(ds, r=0.2, r_max=25, stoyan=0.15, n_rank=1)

# plot PCF
plot(x=pcf, xlim=c(0, 20), ylim=c(0, 2.2))
```


### Links
* [GEOS](https://trac.osgeo.org/geos/)
* [GDAL/OGR Website](https://www.gdal.org/)
* [GDAL/OGR Github Repository](https://github.com/OSGeo/gdal)
* [Rcpp Website](http://www.rcpp.org/)
* [Rcpp Github Repository](https://github.com/RcppCore/Rcpp)
* [R package `sf`, a modern approach to geo data in R](https://github.com/r-spatial/sf)


### References
Nuske, R.S., Sprauer, S. and Saborowski J. (2009): Adapting the pair-correlation function for analysing the spatial distribution of canopy gaps. Forest Ecology and Management 259(1): 107–116. DOI: 10.1016/j.foreco.2009.09.050
