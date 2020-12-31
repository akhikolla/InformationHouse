#' A panel of of 100 plain vanilla fixed coupon corporate bonds.
#'
#' A simulated dataset of 100 plain vanilla fixed coupon
#' corporate bonds issued in 2016.
#'
#' @docType data
#'
#' @usage data(PanelSomeBonds2016)
#'
#' @format A data frame with 12718 rows and 16 variables:
#' \describe{
#'   \item{ID.No}{Identification number of the security.}
#'   \item{Coup.Type}{Type of the bond's coupon.}
#'   \item{Issue.Date}{The bond's issue date. Object of class Date
#'                     with format \code{"\%Y-\%m-\%d"}.}
#'   \item{FIAD.Input}{Date on which the interest accrual starts (so-called
#'               "dated date"). Object of class Date with format
#'               \code{"\%Y-\%m-\%d"}.}
#'   \item{FIPD.Input}{First interest payment date after \code{Issue.Date}.
#'               Object of class Date with format \code{"\%Y-\%m-\%d"}.}
#'   \item{LIPD.Input}{Last interest payment date before \code{Mat.Date}.
#'               Object of class Date with format \code{"\%Y-\%m-\%d"}.}
#'   \item{Mat.Date}{So-called "maturity date" i.e. date on which the
#'                   redemption value and the final interest are paid.
#'                   Object of class Date with format \code{"\%Y-\%m-\%d"}.}
#'   \item{CpY.Input}{Number of interest payments per year. Object of class numeric.}
#'   \item{Coup.Input}{The nominal interest p.a. of the bond in percent. Object
#'               of class numeric.}
#'   \item{RV.Input}{The face value (= redemption value, par value) of
#'             the bond in percent.}
#'   \item{DCC.Input}{The day count convention the bond follows. Type ?AccrInt for details.}
#'   \item{EOM.Input}{Boolean indicating whether the bond follows the End-of-Month rule.}
#'   \item{TradeDate}{The calendar date on which the clean price was observed.}
#'   \item{SETT}{The settlement date that corresponds to \code{TradeDate}.}
#'   \item{CP.Input}{The clean price of the bond on \code{TradeDate}.}
#'   \item{YtM.Input}{The annualized yield to maturity of the bond on \code{TradeDate}.}
#' }
#'
#' @keywords datasets
#'
"PanelSomeBonds2016"
