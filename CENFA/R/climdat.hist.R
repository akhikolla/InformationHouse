#' Historical climate data
#'
#' Historical climate dataset of 10 bioclimate variables (\url{https://www.worldclim.org}).
#'
#' @format A RasterBrick with 10 layers:
#' \describe{
#   \item{MAT}{mean annual temperature (degrees C * 10)}
#'   \item{MDR}{mean diurnal range (mean of monthly max temp - min temp)}
#'   \item{ISO}{isothermality (MDR/TAR * 100)}
#'   \item{TS}{temperature seasonality (sd monthly temp * 100)}
#'   \item{HMmax}{max temp of warmest month}
#'   \item{CMmin}{min temp of coldest month}
#   \item{TAR}{temp annual range (HMmax - CMmin)}
#   \item{MTWQ}{mean temp of wettest quarter}
#   \item{MTDQ}{mean temp of driest quarter}
#   \item{MTHQ}{mean temp of hottest quarter}
#   \item{MTCQ}{mean temp of coldest quarter}
#   \item{TAP}{total annual precipitation (mm)}
#'   \item{PWM}{precip of wettest month}
#'   \item{PDM}{precip of driest month}
#'   \item{PS}{precip seasonality (sd/mean monthly precip)}
#'   \item{PWQ}{precip of wettest quarter}
#'   \item{PDQ}{precip of driest quarter}
#   \item{PHQ}{precip of hottest quarter}
#   \item{PCQ}{precip of coldest quarter}
#' }
#' @source \url{https://www.worldclim.org}
#' @seealso \code{\link{climdat.fut}}, \code{\link{tree-data}}
"climdat.hist"
