# News for Package 'apcf'

## Changes in version 0.1.5
* cleaned up configure
* new maintainer email address due to problem with email provider

## Changes in version 0.1.4
* changed configure to cater to GDAL version 2 and 3 (thanks to package sf)

## Changes in version 0.1.3
* tweaked makevars.win in expectation of changes in win toolchain (thanks @jeroen)

## Changes in version 0.1.2
* added configure.ac to package source bundle
* clarified installation of stable/development version in README
* special treatment for R-devel in configure regarding R version check

## Changes in version 0.1.1
* started using continuous integration tools (Travis, AppVeyor, Drone)
* made it build on windows
* cleanup and wording

## Changes in version 0.1.0
* Intial version. Reimplementation of the Adapted Pair Correlation Function
  in C++ using GEOS and GDAL libraries directly instead of through PostGIS.
  Contains mainly the functions `pat2dists()`, `dists2pcf()`, and 
  a plot-function.
