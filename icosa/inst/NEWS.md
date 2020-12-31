# Change log

## [0.10.0 (build 1016)] - CRAN submission - 2020-02-15
### Changed
- vignette was cut for performance optimization
- documentation was systematically revised


## [0.9.88 (build 1015)] - 2020-02-09
### Added
- documentation for the 'tessguide' dataset
- 'legend' argument to faces3d method of facelayer

### Changed
- the locate() function's main arguments were renamed to 'x' and 'y' with S4 dispatch for 'y'
- the surfacecentroid() function is rewritten with S4 method dispatch, with main argument of 'x'
- corrected roxygen tags for S3 methods
- minor documentation corrections, package ready for CRAN resubmission

## [0.9.88 (build 1014)] - 2020-02-09
### Added
- long at lat arguments to the CarToPol() function's data.frame-method

### Changed
- CarToPol() and PolToCar()  functions rewritten with S4 method dispatch
- main function argument of CarToPol() and PolToCar() is now 'x' instead of 'longLatMat' and 'matXYZ'

## [0.9.88 (build 1013)] - 2020-02-04
### Added
- proper usage (roxygen tags before "function") for triggering help files in case of namespace conflicts.

### Changed
- Reorganization of Roxygen tags: removed all aliases, -method bullshit


## [0.9.88 (build 1012)] - 2020-02-04
### Changed
- rgl package moved to Suggests from Depends, freeing the package from its grasp
- plot3d() rewritten as S3 from S4
- rgl dependent functions have warnings to instruct installation.


## [0.9.88 (build 1011)] - 2020-02-04
### Added
- Collate field for DESCRIPTION

### Changed
- R file structure

## [0.9.88 (build 1010)] - 2020-02-03
### Added
- conditional generic declaration to values<- 
- conditional generic declaration to resample()
- conditional generic declaration to values()
- conditional generic declaration to rotate()

### Changed
- complete NAMESPACE cleanup with roxygen
- package 'raster' moved from 'Imports' to 'Suggests'

### Removed 
- old vignette

----------------------------------------------------------------------


## [0.9.87 (build 1009)] - 2019-08-20
### Added 
- Dynamic sp resolution setting 
- proper NEWS file
- suppressed warnings when triangles are produced
- 'breaks', 'inclusive', 'discrete' arguments added for facelayer plotting method. 

### Changed
- The default colour value of the facelayer plotting function is changed to 'heat'.
- The heatmap generation of the facelayer plotting methods is completely rewritten.
- The default of 'alpha' of the facelayer-plotting argument is changed to NULL. 

### Removed
- the minVal and maxVal argument of heatMapLegend()

### Fixed
- Bug in the heatMapLegend() function that cause the legend to miss the last colour.

## [0.9.86 (build 1008)] - 2019-05-02
### Fixed 
- Proper projection treatment during the 3d plotting of sp-type objects. (Thanks to Dominik Jaskierniak for reporting the error).

### Added
- The 'radius' argument is added to the rgl-extension 3d plotting functions.

### Removed
- The 'inner' argument of the surfacecentroid() and chullsphere() functions is deprecated.

## [0.9.85 (build 1001)] - 2019-03-05
### Added
- the cellocator() function
- added the namedorder argument of the vicinity() function


## [0.9.84] - 2018-12-03
### Added
- the strict argument to the chulshere() function - unfinished!
- the tessellation guide (tessguide) object is added to the package

## [0.9.83] - 2018-10-27
### Added
- 'legend' argument for the facelayer method of 'plot()', to disable the plotting of the legend
- the surfacechullsphere() function to calculate areas of spherical convex hulls

## [0.9.82] - 2018-09-27
Rebuilt from previous version.

### Changed
- vignette structure update
- vignette now uses z1 OpenStreetMap land polygons instead of z3
- the package no longer requires the 'rgdal' package, which is only used for projection changes

### Fixed
- added compatibility for up-to-date Rcpp versions
- explicit self-assignment compilation warning with clang 7.0.0

### Known issues
- warnings produced by sp::Polygons function - will be corrected shortly


## [0.9.81] - 2017-04-18
### Fixed
- memory deallocation issues
- SpPolygons(): zenith/nadir face issue
- unnecessary 'rgdal' namespace import

## [0.9.80] (first Beta) - 2017-04-17
### Added
- the locate() function was updated to version 6.0, now incorporating the 'randomborders' argument
- the occupied() function was rebuilt on the locate() function
- the centers() shorthand was added to the quick extraction of the facecenters
- major work on the help files
- the value replacement method of the gridlayer was extended to host lat/long indicators
- added logical subsetting to the facelayer
- bugfix for the x86 application of SpPolygons()
 
## [0.8.61] - 2017-04-09
### Added
- igraph is added as related content 
- the function gridgraph() is implmented for 'trigrid' and for 'facelayer' classes
- fixed bug in OccupiedFaces() for SpatialPolygons
- the @graph slot has been added to the 'trigrid' class to host 'igraph' class graph representations
- package namespace is cleaned up to include only importing namespaces, only the rgl package is a dependency
- the newgraph() function is added to the package
- added an Rcpp function to interpolate 3d lines
- the gridgraph() function is added for the 'hexagrid' class, and the class constructor is upgraded to include it
- the vicinity() function replaced the neighbours() function
- the class of the values is now printed to the console, when the show method of the facelayer is called
- the tesselation and class of the source grid is now added to the gridlayer class, along with a function that checks the linked grid's compatibility to
the facelayer
- the translate() function is added to the package for fast reposition of translated grids to the origin
- the chullsphere() and surfacecentroid() functions were added to the packages
- resampling methods were added for the 'hexagrid'
- igraph representation was added to the vignette


## [0.8.16] - 2017-02-17
### Added
- the shapes() function is added to the package. It will calculate a value that is proportional to the irregularities of the triangular faces or subfaces.
- character values in facelayers will be plotted with random colors
- rgdal dependecy of spTransorm() was properly resolved
- the missing belts slot of the hexagrid class was filled 
- fixed missing invalid input message for PolToCar() function.
- grid resolution is now displayed with the show() method.
- the group generics Ops, Math, and Summary were implemented for the facelayer.
- basic functions of the 'stats' packages were implemented
- latitude and longitude-wise selection was implemented to the facelayer subsetting method
- resampling of 'facelayer' objects to based on 'trigrid' classes was applied in the program: downscaling
- added 2d plotting scheme for 'facelayers', where the values of the facelayer are exclusively colours

### Deleted
- the deprecated argument "border" is no longer available for the user interface of locate()

## [0.8.0] (Alpha) - 2016-11-25
### Added
- Most gridding features are functional and should be usable, I consider the package to be ready for alpha testing.

### Notes
Pre-alpha versions were not registered.
