# globalVariables(c(".", "!!"))

#' Climate attributes from projections.
#'
#' A example dataset containing the climate attribute values 
#' in fraction/additive change
#'
#' @format A data frame with 6 rows and 6 variables:
#' \describe{
#'   \item{P_ann_tot_m}{change in mean annual total P, fraction}
#'   \item{P_ann_seasRatio_m}{change in seasonal ratio of P, fraction}
#'   \item{P_ann_nWet_m}{change in the number of wet days, fraction}
#'   \item{Temp_ann_avg_m}{change in average annual Temp, additive}
#'   \item{Name}{name of the climate model}
#'   \item{Avg. Deficit}{performance metric values}
#' }
"egClimData"

#' Performance metrics of the tank model using simple scaled scenarios.
#'
#' @format A list with 2 elements
#' \describe{
#'   \item{Avg. Deficit}{average daily deficit of water, litres}
#'   \item{Reliability}{reliability of the tank, fraction}
#' }
"egScalPerformance"

#' Summary of a simple scaled scenario.
#'
#' Summary generated using the function \code{getSimSummary}.
#'
#' @format A list containing 3 elements
#' \describe{
#'   \item{simDates}{the dates of the simulation}
#'   \item{expSpace}{the exposure space of the simulation}
#'   \item{controlFile}{"scaling"}
#' }
"egScalSummary"

#' Summary of a regGrid scenario.
#'
#' Summary generated using the function \code{getSimSummary} for
#' a scenarios generated using stochastic models for a regGrid exposure space
#'
#' @format A list containing 13 elements
"egSimSummary"


#' Summary of a OAT scenario.
#'
#' Summary generated using the function \code{getSimSummary} for
#' a scenarios generated using stochastic models for an OAT exposure space
#'
#' @format A list containing 13 elements
"egSimOATSummary"

#' Performance metrics of the tank model using OAT scenarios.
#'
#' @format A list with 2 elements
#' \describe{
#'   \item{Avg. Deficit}{average daily deficit of water, litres}
#'   \item{Reliability}{reliability of the tank, fraction}
#' }
"egSimOATPerformance"

#' Performance metrics of the tank model using regGrid scenarios.
#'
#' @format A list with 2 elements
#' \describe{
#'   \item{Avg. Deficit}{average daily deficit of water, litres}
#'   \item{Reliability}{reliability of the tank, fraction}
#' }
"egSimPerformance"

#' Performance metrics of an alternate tank model using regGrid scenarios.
#'
#' @format A list with 2 elements
#' \describe{
#'   \item{Avg. Deficit}{average daily deficit of water, litres}
#'   \item{Reliability}{reliability of the tank, fraction}
#' }
"egSimPerformance_systemB"


