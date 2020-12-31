#' Example data for calculating bouts of moderate-to-vigorous physical activity
#'
#' A dataset containing accelerometer data and predicted energy expenditure in
#' metabolic equivalents (METs) that can be used to classify
#' moderate-to-vigorous physical activity in continuous bouts.
#'
#' @format A data frame with 10080 rows and 12 variables:
#' \describe{
#'   \item{FileID}{character. Name of the file originating the data}
#'   \item{Date}{character giving the date ("\%m/\%d/\%Y")}
#'   \item{Time}{character giving the time ("\%H:\%M:\%S")}
#'   \item{DateTime}{full timestamp (\%Y-\%m-\%d \%H:\%M:\%S) given as character}
#'   \item{dayofyear}{numeric day of the year (i.e., julian date)}
#'   \item{minofday}{numeric minute of the day (i.e., 0 for midnight and 1439
#'   for 11:59)}
#'   \item{Axis1}{activity counts for the device's first axis}
#'   \item{Axis2}{activity counts for the device's second axis}
#'   \item{Axis3}{activity counts for the device's third axis}
#'   \item{Steps}{number of steps taken}
#'   \item{Vector.Magnitude}{vector magnitude (Euclidian norm) of the activity
#'   counts from the three axes}
#'   \item{METs}{predicted energy expenditure, in metabolic equivalents}
#' }
"ex_data"
