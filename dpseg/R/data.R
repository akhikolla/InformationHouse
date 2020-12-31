### DATA SET OD growth curves by Tom Rohr, Anna Behle, 2018

#' Escherichia coli growth curves.
#'
#' Optical density (OD) data from a 96-well microtiter plate experiment,
#' growing Escherichia coli cells in M9 medium in a BMG Optima platereader.
#' 
#' @source Tom Rohr, Anna Behle, Rainer Machne, HHU Duesseldorf, 2018
#' @format A data frame with the measurement time in column 1 and
#'     bacterial growth data (or blanks) in
#'     \code{2:ncol(oddata)}. Column names correspond to the well on
#'     the microtiter plate.
"oddata"

### RNA-seq read-count data
### 
### @source Douglas B. Murray, Rainer Machne
### @format A data frame with columns:
### \describe{
### \item{position:}{genomic position}
### \item{readcount:}{RNA-seq read-counts}
### }
 ##"rnaseq"
NULL
