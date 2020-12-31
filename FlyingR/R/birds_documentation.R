#' Sample 28 birds
#'
#' Preset birds data, extracted from Flight Pennycuick(2008). Fat mass percentage
#' generated randomly where zero.
#'
#' @format A data frame with 28 observations and 5 variables not counting the
#'         name.
#'
#' \describe{
#'    \item{Scientific.name}{Name of bird species}
#'    \item{Empty.mass}{Body mass in Kg. Includes fuel (fat mass). In this case the crops
#'    were empty but otherwise one should always use the all-up mass (body mass + crop)}
#'    \item{Wing.span}{Length of wings spread out in meters}
#'    \item{Fat.mass}{Mass of fat that is consumable as fuel in Kg}
#'    \item{Order}{Order of the species (passerine vs non-passerine)}
#'    \item{Wing.area}{Area of both wing projected on a flat surface in meters
#'    squared}
#'    \item{Muscle.mass}{Mass in Kg. of flight muscles}
#'
#' }
#'
"birds"
