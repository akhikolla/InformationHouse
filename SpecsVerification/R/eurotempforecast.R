#' Seasonal ensemble forecast of European average summer temperature
#'
#' A hindcast dataset of average European (30N,75N,12.5W,42.5E) summer (June/July/August) surface temperatures. Forecasts were initialised in May the same year. Observations and 15-member ensemble forecasts were derived from the publicly available NCEP Reanalysis (Suranjana, 2010) and the NCEP Climate Forecast System Version 2 (Suranjana, 2014), respectively. The data was downloaded through the ECOMS User Data Gateway (Santander Meteorology Group, 2015).
#'
#' @format Variables contained in the data set:
#' \itemize{
#'   \item `obs` average European summer temperature observations 
#'   \item `ens` mean-debiased ensemble forecast data, i.e. mean(ens) == mean(obs)
#'   \item `obs.lag` the observations lagged by one year, same length as `obs`
#'   \item `obs.bin` binary observations (0 or 1), obs[i] = 1 indicates that the temperature of year i exceeded the temperature of year i-1
#'   \item `ens.bin` binary ensemble forecast (each member is either 0 or 1), ens[i, j] = 1 if the j-th ensemble member in year i exceeded the observed temperature of year i-1
#'   \item `obs.cat` categorical observations. obs.cat[i] is either 1, 2, and 3, indicating that the temperature in year i was lower, similar, higher than temperature in year i-1. Similar is defined as within a half degree interval centered around last years temperature.
#'   \item `ens.cat` categorical ensemble forecast. ens.cat[i, j] is either 1, 2, or 3. The categories are defined as for `obs.cat`.
#' }
#' @docType data
#' @keywords datasets
#' @name eurotempforecast
#' @aliases ens obs obs.lag ens.bin obs.bin ens.cat obs.cat
#' @references 
#' Saha, Suranjana, and Coauthors, 2010: The NCEP Climate Forecast System Reanalysis. Bull. Amer. Meteor. Soc., 91, 1015.1057. \doi{10.1175/2010BAMS3001.1}
#' Saha, Suranjana and Coauthors, 2014: The NCEP Climate Forecast System Version 2. J. Clim., 27, 2185--2208, \doi{10.1175/JCLI-D-12-00823.1}
#' Santander Meteorology Group (2015). ecomsUDG.Raccess: R interface to the ECOMS User Data Gateway. R package version 4.2-0. \url{http://meteo.unican.es/trac/wiki/udg/ecoms}
#' @usage data(eurotempforecast)
NULL

