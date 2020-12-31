#' AnnivDates (time-invariant properties and temporal structure)
#'
#' \bold{AnnivDates} returns a bond's time-invariant characteristics and temporal structure as a list of
#' three or four named data frames.
#'
#' \bold{AnnivDates} generates a list of the three data frames \code{Warnings}, \code{Traits}
#' and \code{DateVectors}. If the variable \code{Coup} is passed to the function,
#' the output contains additionally the data frame \code{PaySched}. \bold{AnnivDates} is meant to analyze
#' large data frames. Therefore some features are implemented to evaluate the quality of the data. The
#' output of these features is stored in the data frame \code{Warnings}. Please see section \bold{Value}
#' for a detailed description of the tests run and the meaning of the variables in \code{Warnings}. The
#' data frame \code{Traits} contains all time-invariant bond characteristics that were either provided by
#' the user or calculated by the function. The data frame \code{DateVectors} contains three vectors
#' of Date-Objects named \code{RealDates}, \code{CoupDates} and \code{AnnivDates} and three vectors of
#' numerics named \code{RD_indexes}, \code{CD_indexes} and \code{AD_indexes}. These vectors are
#' used in the other functions of this package according to the methodology presented in Djatschenko (2018).
#' The data frame \code{PaySched} matches \code{CoupDates}
#' to the actual amount of interest that the bond pays on the respective interest payment date. Section
#' \bold{Value} provides further information on the output of the function \bold{AnnivDates}. Below
#' information on the proper input format is provided. Subsequently follows information on the operating
#' principle of the function \bold{AnnivDates} and on the assumptions that are met to
#' estimate the points in time needed to evaluate a bond.
#'
#' \itemize{
#'   \item The dates \code{Em}, \code{Mat}, \code{FIPD}, \code{LIPD} and \code{FIAD} can be provided as
#'   \enumerate{
#'   \item "Date" with format \code{"\%Y-\%m-\%d"}, or
#'   \item "numeric" with the appropriate \code{DateOrigin}, or
#'   \item number of class "character" with the appropriate \code{DateOrigin}, or
#'   \item string of class "character" in the format \code{"yyyy-mm-dd"}.
#'   }
#'   \code{CpY}, \code{RV} and \code{Coup} can be provided either as class "numeric" or as a number of
#'   class "character".
#'
#'   \item The provided issue date (\code{Em}) is instantly substituted by the first interest accrual
#'   date (\code{FIAD}) if \code{FIAD} is available and different from \code{Em}.
#'
#'   \item Before the determination of the bond's date characteristics begins, the code evaluates
#'   the provided calendar dates for plausibility. In this process implausible dates are dropped.
#'   The sort of corresponding implausibility is identified and stored in a warning flag. (See
#'   section \bold{Value} for details.)
#'
#'   \item The remaining valid calendar dates are used to gauge whether the bond follows the
#'   End-of-Month-Rule. The resulting parameter est_EOM can take on the following values:
#'      \describe{
#'        \item{}{
#'          \tabular{cl}{
#'          \bold{\emph{Case 1:}} \tab \bold{\code{FIPD} and \code{LIPD} are both \code{NA}} \cr
#'                    ___________ \tab ____________________________________                  \cr
#'             \code{est_EOM = 1} \tab , if \code{Mat} is the last day of a month.           \cr
#'             \code{est_EOM = 0} \tab , else.                                               \cr
#'                     ========== \tab ================================                      \cr
#'          }
#'        }
#'        \item{}{
#'          \tabular{cl}{
#'          \bold{\emph{Case 2:}} \tab \bold{\code{FIPD} is \code{NA} and \code{LIPD} is a valid calendar date} \cr
#'                    ___________ \tab ____________________________________                                     \cr
#'             \code{est_EOM = 1} \tab , if \code{LIPD} is the last day of a month.                             \cr
#'             \code{est_EOM = 0} \tab , else.                                                                  \cr
#'                     ========== \tab ================================                                         \cr
#'          }
#'        }
#'        \item{}{
#'          \tabular{cl}{
#'          \bold{\emph{Case 3:}} \tab \bold{\code{FIPD} is a valid calendar date and \code{LIPD} is \code{NA}} \cr
#'                    ___________ \tab ____________________________________                                     \cr
#'             \code{est_EOM = 1} \tab , if \code{FIPD} is the last day of a month.                             \cr
#'             \code{est_EOM = 0} \tab , else.                                                                  \cr
#'                     ========== \tab ================================                                         \cr
#'          }
#'        }
#'        \item{}{
#'          \tabular{cl}{
#'          \bold{\emph{Case 4:}} \tab \bold{\code{FIPD} and \code{LIPD} are valid calendar dates} \cr
#'                    ___________ \tab ____________________________________                        \cr
#'             \code{est_EOM = 1} \tab , if \code{LIPD} is the last day of a month.                \cr
#'             \code{est_EOM = 0} \tab , else.                                                     \cr
#'                     ========== \tab ================================                            \cr
#'          }
#'        }
#'      }
#'   \item If \code{EOM} is initially missing or \code{NA} or not element of \code{\{0,1\}}, \code{EOM}
#'   is set \code{est_EOM} with a warning.
#'   \item If the initially provided value of \code{EOM} deviates from \code{est_EOM}, the following two
#'   cases apply:
#'      \tabular{cl}{
#'              ________ \tab _________________________________________                       \cr
#'               Case 1: \tab If \code{EOM = 0} and \code{est_EOM = 1}:                       \cr
#'                       \tab \code{EOM} is not overridden and remains \code{EOM = 0}         \cr
#'              ________ \tab _________________________________________                       \cr
#'               Case 2: \tab If \code{EOM = 1} and \code{est_EOM = 0}:                       \cr
#'                       \tab \code{EOM} is overridden and set \code{EOM = 0} with a warning. \cr
#'                       \tab Keeping \code{EOM = 1} in this case would conflict with         \cr
#'                       \tab the provided \code{Mat}, \code{FIPD} or \code{LIPD}.            \cr
#'              ________ \tab _________________________________________                       \cr
#'                 Note: \tab Set the option \code{FindEOM=TRUE} to always use                \cr
#'                       \tab \code{est_EOM} found by the code.                               \cr
#'               ======= \tab ====================================                            \cr
#'               }
#'
#'   \item If \code{FIPD} and \code{LIPD} are both available, the lengths of the first and final coupon
#'   periods are determinate and can be "regular", "long" or "short". To find the interest payment dates
#'   between \code{FIPD} and \code{LIPD} the following assumptions are met:
#'   \enumerate{
#'   \item \preformatted{The interest payment dates between FIPD and LIPD are
#'   evenly distributed.}
#'   \item \preformatted{The value of EOM determines the location of
#'   all interest payment dates.}
#'   }
#'   If assumption 1 is violated, the exact locatations of the interest payment dates between
#'   \code{FIPD} and \code{LIPD} are ambiguous. The assumption is violated particularly, if
#'   \enumerate{
#'   \item \code{FIPD} and \code{LIPD} are in the same month of the same year but not on the same day, or
#'   \item the month difference between \code{FIPD} and \code{LIPD} is not a multiple of the number
#'   of months implied by \code{CpY}, or
#'   \item \code{FIPD} and \code{LIPD} are not both last day in month,
#'         their day figures differ and the day figure difference between \code{FIPD}
#'         and \code{LIPD} is not due to different month lengths.
#'   }
#'   In each of the three cases, \code{FIPD} and \code{LIPD} are dropped
#'   with the flag \code{IPD_CpY_Corrupt = 1}.
#'
#'   \item If neither \code{FIPD} nor \code{LIPD} are available the code
#'   evaluates the bond based only upon the required variables \code{Em} and
#'   \code{Mat} (and \code{CpY}, which is \code{2} by default). Since FIPD is
#'   not given, it is impossible to distinguish between a "short" and "long" odd
#'   first coupon period, without an assumption on the number of interest
#'   payment dates. Consequently the first coupon period is assumed to be either
#'   "regular" or "short". The locations of \code{FIPD} and \code{LIPD} are
#'   estimated under the following assumptions:
#'   \enumerate{
#'   \item \preformatted{The final coupon period is "regular".}
#'   \item \preformatted{The interest payment dates between the estimated
#'   FIPD and Mat are evenly distributed.}
#'   \item \preformatted{The value of EOM determines the location of
#'   all interest payment dates.}
#'   }
#'
#'   \item If \code{LIPD} is available but \code{FIPD} is not, the length
#'   of the final coupon payment period is determined by \code{LIPD} and
#'   \code{Mat} and can be "regular", "long" or "short". The locations of
#'   the interest payment dates are estimated under the following assumptions:
#'   \enumerate{
#'   \item \preformatted{The first coupon period is either "regular" or "short".}
#'   \item \preformatted{The interest payment dates between the estimated
#'   FIPD and LIPD are evenly distributed.}
#'   \item \preformatted{The value of EOM determines the location of
#'   all interest payment dates.}
#'   }
#'
#'   \item If \code{FIPD} is available but \code{LIPD} is not, the length
#'   of the first coupon payment period is determined by \code{Em} and
#'   \code{FIPD} and can be "regular", "long" or "short". The locations of
#'   the interest payment dates are estimated under the following assumptions:
#'   \enumerate{
#'   \item \preformatted{The final coupon period is either "regular" or "short".}
#'   \item \preformatted{The interest payment dates between FIPD and
#'   the estimated LIPD are evenly distributed.}
#'   \item \preformatted{The value of EOM determines the location of
#'   all interest payment dates.}
#'   }
#' }
#'
#'
#' @param Em The bond's issue date. (required)
#' @param Mat Maturity date, i.e. date on which the redemption value and the final interest
#'        are paid. (required)
#' @param CpY Number of interest payments per year (non-negative integer; element of the set
#'        \{0,1,2,3,4,6,12\}. Default: 2.
#' @param FIPD First interest payment date after \code{Em}.
#' @param LIPD Last interest payment date prior to \code{Mat}.
#' @param FIAD Date on which the interest accrual starts (so-called "dated date").
#' @param RV The redemption value of the bond. Default: 100.
#' @param Coup Nominal interest rate per year in percent. Default: \code{NA}.
#' @param DCC The day count convention the bond follows. Default: \code{NA}.
#'        For a list of day count conventions currently implemented type \code{View(List.DCC)}.
#' @param EOM Boolean indicating whether the bond follows the End-of-Month rule. Default: \code{NA}.
#' @param DateOrigin Determines the starting point for the daycount in "Date" objects.
#'        Default: "1970-01-01".
#' @param InputCheck If 1, the input variables are checked for the correct format. Default: 1.
#' @param FindEOM If \code{TRUE}, \code{EOM} is overridden by the value inferred from the data.
#'        Default: \code{FALSE}.
#' @param RegCF.equal If 0, the amounts of regular cash flows are calculated according to the
#'        stipulated \code{DCC}. Any other value forces all regular cash flows to be equal sized.
#'        Default: 0.
#'
#' @return All dates are returned irrespective of whether they are on a business day or not.
#'  \describe{
#'    \item{\emph{\bold{DateVectors}} (data frame)}{
#'      \describe{
#'        \item{}{}
#'        \item{\emph{RealDates}}{A vector of Date class objects with format "\%Y-\%m-\%d" in ascending order,
#'        that contains the issue date, all actual coupon payment dates and the maturity date.}
#'        \item{\emph{RD_indexes}}{A vector of numerics capturing the temporal structure of the bond.}
#'        \item{\emph{CoupDates}}{A vector of Date class objects with format "\%Y-\%m-\%d" in ascending order,
#'        that contains all actual coupon payment dates and the maturity date.}
#'        \item{\emph{CD_indexes}}{A vector of numerics capturing the temporal structure of the bond.}
#'        \item{\emph{AnnivDates}}{A vector of Date class objects with format "\%Y-\%m-\%d" in ascending order,
#'        that contains all theoretical coupon anniversary dates. The first value of \emph{AnnivDates} is the
#'        anniversary date immediately preceding the issue date, if the bond has an irregular first coupon
#'        period; otherwise it is the issue date. The final value of \emph{AnnivDates} is the anniversary
#'        date immediately succeeding the maturity date, if the bond has an irregular final coupon period;
#'        otherwise it is the maturity date.}
#'        \item{\emph{AD_indexes}}{A vector of numerics capturing the temporal structure of the bond.}
#'        \item{}{}
#'      }
#'    }
#'    \item{\emph{\bold{PaySched}} (data frame)}{
#'      \describe{
#'        \item{}{}
#'        \item{\emph{CoupDates}}{A vector of Date class objects with format "\%Y-\%m-\%d" in ascending order,
#'        that contains all actual coupon payment dates and the maturity date.}
#'        \item{\emph{CoupPayments}}{A vector of class "numeric" objects, that contains the actual amounts of
#'        interest that the bond pays on the respective coupon payment dates. The unit of these payments is the
#'        same as that of \code{RV} that was passed to the function. \code{RV} is not included in the final
#'        interest payment.}
#'        \item{\bold{NOTE:}}{\code{PaySched} is created only if the variable \code{Coup} is provided.}
#'        \item{}{}
#'      }
#'    }
#'    \item{\emph{\bold{Traits}} (data frame)}{
#'      \describe{
#'        \item{}{}
#'        \item{\emph{DateOrigin}}{The starting point for the daycount in "Date" objects.}
#'        \item{\emph{CpY}}{Number of interest payments per year.}
#'        \item{\emph{FIAD}}{Date on which the interest accrual starts (so-called "dated date").}
#'        \item{\emph{Em}}{The bond's issue date that was used for calculations.}
#'        \item{\emph{Em_Orig}}{The bond's issue date that was entered.}
#'        \item{\emph{FIPD}}{The first interest payment date after \code{Em} that was used for calculations.
#'                           If the entered \code{FIPD} was dropped during the calculation process,
#'                           the value is \code{NA}.}
#'        \item{\emph{FIPD_Orig}}{The first interest payment date after \code{Em} that was entered.}
#'        \item{\emph{est_FIPD}}{The estimated first interest payment date after \code{Em}. \code{NA}, if
#'                               a valid \code{FIPD} was entered.}
#'        \item{\emph{LIPD}}{The last interest payment date prior to \code{Mat} that was used for
#'                           calculations. If the entered \code{LIPD} was dropped during the calculation
#'                          process, the value is \code{NA}.}
#'        \item{\emph{LIPD_Orig}}{The last interest payment date prior to \code{Mat} that was entered.}
#'        \item{\emph{est_LIPD}}{The estimated last interest payment date prior to \code{Mat}. \code{NA},
#'                               if a valid \code{LIPD} was entered.}
#'        \item{\emph{Mat}}{The maturity date that was entered.}
#'        \item{\emph{Refer}}{Reference date that determines the day figures of all AnnivDates.}
#'        \item{\emph{FCPType}}{A character string indicating the type of the first coupon period.
#'        Values: "long", "regular", "short".}
#'        \item{\emph{FCPLength}}{Length of the first coupon period as a fraction of a regular coupon period.}
#'        \item{\emph{LCPType}}{A character string indicating the type of the last coupon period.
#'        Values: "long", "regular", "short".}
#'        \item{\emph{LCPLength}}{Length of the final coupon period as a fraction of a regular coupon period.}
#'        \item{\emph{Par}}{The redemption value of the bond.}
#'        \item{\emph{CouponInPercent.p.a}}{Nominal interest rate per year in percent.}
#'        \item{\emph{DayCountConvention}}{The day count convention the bond follows.}
#'        \item{\emph{EOM_Orig}}{The value of \code{EOM} that was entered.}
#'        \item{\emph{est_EOM}}{The estimated value of \code{EOM}.}
#'        \item{\emph{EOM_used}}{The value of \code{EOM} that was used in the calculations.}
#'        \item{}{}
#'      }
#'    }
#'    \item{\emph{\bold{Warnings}} (data frame)}{
#'      \describe{
#'        \item{}{}
#'        \item{}{A set of flags that indicate the occurrence of warnings during the execution.
#'                Below they are listed according to the hierarchical structure within the function \bold{AnnivDates}.}
#'        \item{}{
#'          \tabular{rcl}{
#'                                            \tab     \tab                                                                         \cr
#'                                            \tab     \tab                                                                         \cr
#'                                            \tab     \tab                                                                         \cr
#'            \bold{\emph{Em_FIAD_differ}} =  \tab     \tab                                                                         \cr
#'                                          1 \tab     \tab , if the provided issue date (\code{Em}) was substituted by the first   \cr
#'                                            \tab     \tab interest accrual date (\code{FIAD}).                                    \cr
#'                                            \tab     \tab This happens, if \code{FIAD} is available and different from \code{Em}. \cr
#'                                            \tab     \tab ________________________________________________                        \cr
#'                                            \tab     \tab \emph{Note:} No warning is displayed.                                   \cr
#'                        ___________________ \tab ___ \tab ________________________________________________                        \cr
#'                                          0 \tab     \tab , else.                                                                 \cr
#'                          ================= \tab === \tab ===========================================                             \cr
#'          }
#'        }
#'        \item{}{
#'          \tabular{rcl}{
#'            \bold{\emph{EmMatMissing}} =  \tab     \tab                                                                                        \cr
#'                                        1 \tab     \tab , if either issue date (\code{Em}) or maturity date (\code{Mat}) or both               \cr
#'                                          \tab     \tab are missing or \code{NA}.                                                              \cr
#'                                          \tab     \tab ________________________________________________                                       \cr
#'                                          \tab     \tab \emph{Output:} \emph{RealDates} \code{= NA}, \emph{CoupDates} \code{= NA},             \cr
#'                                          \tab     \tab \emph{AnnivDates} \code{= NA}, \emph{FCPType} \code{= NA}, \emph{LCPType} \code{= NA}. \cr
#'                      ___________________ \tab ___ \tab ________________________________________________                                       \cr
#'                                        0 \tab     \tab , else.                                                                                \cr
#'                        ================= \tab === \tab ===========================================                                            \cr
#'          }
#'        }
#'        \item{}{
#'          \tabular{rcl}{
#'            \bold{\emph{CpYOverride}} =  \tab     \tab                                                                                \cr
#'                                       1 \tab     \tab , if number of interest periods per year (\code{CpY}) is missing or            \cr
#'                                         \tab     \tab \code{NA}, or if the provided \code{CpY} is not element of \{0,1,2,3,4,6,12\}. \cr
#'                                         \tab     \tab ________________________________________________                               \cr
#'                                         \tab     \tab \emph{Note:} \code{CpY} is set 2, and the execution continues.                 \cr
#'                                         \tab     \tab ________________________________________________                               \cr
#'                                         \tab     \tab \emph{Output:} as if \code{CpY} = 2 was provided initially.                    \cr
#'                     ___________________ \tab ___ \tab ________________________________________________                               \cr
#'                                       0 \tab     \tab , else.                                                                        \cr
#'                       ================= \tab === \tab ===========================================                                    \cr
#'          }
#'        }
#'        \item{}{
#'          \tabular{rcl}{
#'            \bold{\emph{RV_set100percent}} =  \tab     \tab                                                                 \cr
#'                                            1 \tab     \tab , if the redemption value (\code{RV}) is missing or \code{NA}.  \cr
#'                                              \tab     \tab ________________________________________________                \cr
#'                                              \tab     \tab \emph{Note:} \code{RV} is set 100, and the execution continues. \cr
#'                                              \tab     \tab ________________________________________________                \cr
#'                                              \tab     \tab \emph{Output:} as if \code{RV} = 100 was provided initially.    \cr
#'                          ___________________ \tab ___ \tab ________________________________________________                \cr
#'                                            0 \tab     \tab , else.                                                         \cr
#'                            ================= \tab === \tab ===========================================                     \cr
#'          }
#'        }
#'        \item{}{
#'          \tabular{rcl}{
#'            \bold{\emph{NegLifeFlag}} =  \tab     \tab                                                                                        \cr
#'                                       1 \tab     \tab , if the provided maturity date (\code{Mat}) is before or on the                       \cr
#'                                         \tab     \tab provided issue date (\code{Em}).                                                       \cr
#'                                         \tab     \tab ________________________________________________                                       \cr
#'                                         \tab     \tab \emph{Output:} \emph{RealDates} \code{= NA}, \emph{CoupDates} \code{= NA},             \cr
#'                                         \tab     \tab \emph{AnnivDates} \code{= NA}, \emph{FCPType} \code{= NA}, \emph{LCPType} \code{= NA}. \cr
#'                     ___________________ \tab ___ \tab ________________________________________________                                       \cr
#'                                       0 \tab     \tab , else.                                                                                \cr
#'                       ================= \tab === \tab ===========================================                                            \cr
#'          }
#'        }
#'        \item{}{
#'          \tabular{rcl}{
#'            \bold{\emph{ZeroFlag}} =  \tab     \tab                                                                                              \cr
#'                                    1 \tab     \tab , if number of interest payments per year (\code{CpY}) is \code{0}.                          \cr
#'                                      \tab     \tab ________________________________________________                                             \cr
#'                                      \tab     \tab \emph{Output:} \emph{RealDates} \code{= (Em,Mat)}, \emph{CoupDates} \code{= Mat},            \cr
#'                                      \tab     \tab \emph{AnnivDates} \code{= (Em,Mat)}, \emph{FCPType} \code{= NA}, \emph{LCPType} \code{= NA}. \cr
#'                  ___________________ \tab ___ \tab ________________________________________________                                             \cr
#'                                    0 \tab     \tab , else.                                                                                      \cr
#'                    ================= \tab === \tab ===========================================                                                  \cr
#'          }
#'        }
#'        \item{}{
#'          \tabular{rcl}{
#'            \bold{\emph{Em_Mat_SameMY}} =  \tab     \tab                                                                                   \cr
#'                                         1 \tab     \tab , if the issue date (\code{Em}) and the maturity date (\code{Mat}) are in the     \cr
#'                                           \tab     \tab same month of the same year but not on the same day, while                        \cr
#'                                           \tab     \tab \code{CpY} is an element of \{1,2,3,4,6,12\}.                                     \cr
#'                                           \tab     \tab ________________________________________________                                  \cr
#'                                           \tab     \tab \emph{Output:} \emph{RealDates} \code{= (Em,Mat)}, \emph{CoupDates} \code{= Mat}, \cr
#'                                           \tab     \tab \emph{FCPType} \code{= short}, \emph{LCPType} \code{= short}.                     \cr
#'                       ___________________ \tab ___ \tab ________________________________________________                                  \cr
#'                                         0 \tab     \tab , else.                                                                           \cr
#'                         ================= \tab === \tab ===========================================                                       \cr
#'          }
#'        }
#'        \item{}{
#'          \tabular{rcl}{
#'            \bold{\emph{ChronErrorFlag}} =  \tab     \tab                                                                                     \cr
#'                                          1 \tab     \tab , if the provided dates are in a wrong chronological order.                         \cr
#'                                            \tab     \tab ________________________________________________                                    \cr
#'                                            \tab     \tab \emph{Note:}                                                                        \cr
#'                                            \tab     \tab The correct ascending chronological order is:                                       \cr
#'                                            \tab     \tab issue date (\code{Em}), first interest payment date (\code{FIPD}),                  \cr
#'                                            \tab     \tab last interest payment date (\code{LIPD}), maturity date (\code{Mat}).               \cr
#'                                            \tab     \tab \code{FIPD} and \code{LIPD} are set \code{as.Date(NA)}.                             \cr
#'                                            \tab     \tab ________________________________________________                                    \cr
#'                                            \tab     \tab \emph{Output:} as if \code{FIPD} and \code{LIPD} were not provided initially.       \cr
#'                        ___________________ \tab ___ \tab ________________________________________________                                    \cr
#'                                          0 \tab     \tab , else.                                                                             \cr
#'                          ================= \tab === \tab ===========================================                                         \cr
#'          }
#'        }
#'        \item{}{
#'          \tabular{rcl}{
#'            \bold{\emph{FIPD_LIPD_equal}} =  \tab     \tab                                                                                           \cr
#'                                           1 \tab     \tab if \code{Em} < \code{FIPD} = \code{LIPD} < \code{Mat}.                                    \cr
#'                                             \tab     \tab ________________________________________________                                          \cr
#'                                             \tab     \tab \emph{Output:} \emph{AnnivDates} contains \code{FIPD} and has at least \code{3} elements. \cr
#'                                             \tab     \tab \emph{RealDates} \code{= (Em,FIPD,Mat)}, \emph{CoupDates} \code{= (FIPD,Mat)}.            \cr
#'                                             \tab     \tab \emph{FCPType} and \emph{LCPType} can be "short", "regular" or "long".                    \cr
#'                         ___________________ \tab ___ \tab ________________________________________________                                          \cr
#'                                           0 \tab     \tab , else.                                                                                   \cr
#'                           ================= \tab === \tab ===========================================                                               \cr
#'          }
#'        }
#'        \item{}{
#'          \tabular{rcl}{
#'            \bold{\emph{IPD_CpY_Corrupt}} =  \tab     \tab                                                                       \cr
#'                                           1 \tab     \tab , if the provided first interest payment date (\code{FIPD}) and last  \cr
#'                                             \tab     \tab interest payment date (\code{LIPD}) are inconsistent with the         \cr
#'                                             \tab     \tab provided number of interest payments per year (\code{CpY}).           \cr
#'                                             \tab     \tab ________________________________________________                      \cr
#'                                             \tab     \tab \emph{Note:}                                                          \cr
#'                                             \tab     \tab Inconsistency occurs if                                               \cr
#'                                             \tab     \tab 1. \code{FIPD} and \code{LIPD} are in the same month of the same year \cr
#'                                             \tab     \tab    but not on the same day, or                                        \cr
#'                                             \tab     \tab 2. the number of months between \code{FIPD} and \code{LIPD} is not a  \cr
#'                                             \tab     \tab    multiple of the number of months implied by \code{CpY}, or         \cr
#'                                             \tab     \tab 3. \code{FIPD} and \code{LIPD} are not both last day in month, their  \cr
#'                                             \tab     \tab    day figures differ and the day figure difference between           \cr
#'                                             \tab     \tab    \code{FIPD} and \code{LIPD} is not due to different month lengths. \cr
#'                                             \tab     \tab                                                                       \cr
#'                                             \tab     \tab In each of the three cases keeping the provided values of              \cr
#'                                             \tab     \tab \code{FIPD} and \code{LIPD} would violate the assumption, that the    \cr
#'                                             \tab     \tab anniversary dates between \code{FIPD} and \code{LIPD} are evenly      \cr
#'                                             \tab     \tab distributed.                                                          \cr
#'                                             \tab     \tab ________________________________________________                      \cr
#'                                             \tab     \tab \code{FIPD} and \code{LIPD} are set \code{as.Date(NA)}                \cr
#'                                             \tab     \tab and the execution continues.                                          \cr
#'                                             \tab     \tab ________________________________________________                      \cr
#'                                             \tab     \tab \emph{Output:}                                                        \cr
#'                                             \tab     \tab as if \code{FIPD} and \code{LIPD} were not provided initially.        \cr
#'                         ___________________ \tab ___ \tab ________________________________________________                      \cr
#'                                           0 \tab     \tab , else.                                                               \cr
#'                           ================= \tab === \tab ===========================================                           \cr
#'          }
#'        }
#'        \item{}{
#'          \tabular{rcl}{
#'            \bold{\emph{EOM_Deviation}} =  \tab     \tab                                                                                  \cr
#'                                         1 \tab     \tab , if the provided value of \code{EOM} deviates from the value that               \cr
#'                                           \tab     \tab is inferred from the provided calendar dates.                                    \cr
#'                                           \tab     \tab ________________________________________________                                 \cr
#'                                           \tab     \tab \emph{Note:}                                                                     \cr
#'                                           \tab     \tab The program analyses the valid values of \code{Em}, \code{Mat}, \code{FIPD} and  \cr
#'                                           \tab     \tab \code{LIPD} to determine the appropriate value of \code{EOM}.                    \cr
#'                                           \tab     \tab                                                                                  \cr
#'                                           \tab     \tab If the initially provided value of \code{EOM} deviates from the value            \cr
#'                                           \tab     \tab determined by the program, there might be an inconsistency                       \cr
#'                                           \tab     \tab in the provided data.                                                            \cr
#'                       ___________________ \tab ___ \tab ________________________________________________                                 \cr
#'                                         0 \tab     \tab , else.                                                                          \cr
#'                         ================= \tab === \tab ===========================================                                      \cr
#'                                           \tab     \tab                                                                                  \cr
#'          }
#'        }
#'        \item{}{
#'          \tabular{rcl}{
#'            \bold{\emph{EOMOverride}} =  \tab     \tab                                                                               \cr
#'                                       1 \tab     \tab , if the provided value of \code{EOM} is overridden by a value that           \cr
#'                                         \tab     \tab is inferred from the provided calendar dates.                                 \cr
#'                                         \tab     \tab ________________________________________________                              \cr
#'                                         \tab     \tab \emph{Note:}                                                                  \cr
#'                                         \tab     \tab This happens automatically if \code{EOM} is initially missing or \code{NA}    \cr
#'                                         \tab     \tab or not element of \code{\{0,1\}} and if the provided value of \code{EOM}      \cr
#'                                         \tab     \tab conflicts with the provided values of \code{FIPD}, \code{LIPD} or \code{Mat}, \cr
#'                                         \tab     \tab e.g. if \code{est_EOM = 0} but \code{EOM = 1}.                                \cr
#'                                         \tab     \tab If \code{EOM_Deviation = 1} and the option \code{FindEOM} is set \code{TRUE}, \cr
#'                                         \tab     \tab the initially provided value of \code{EOM} is also overridden by the          \cr
#'                                         \tab     \tab value that is inferred from the provided calendar dates if                    \cr
#'                                         \tab     \tab \code{est_EOM = 1} but \code{EOM = 0}.                                        \cr
#'                                         \tab     \tab ________________________________________________                              \cr
#'                                         \tab     \tab \emph{Output:}                                                                \cr
#'                                         \tab     \tab as if the value of \code{EOM} that is found by the program was                \cr
#'                                         \tab     \tab provided initially.                                                           \cr
#'                     ___________________ \tab ___ \tab ________________________________________________                              \cr
#'                                       0 \tab     \tab , else.                                                                       \cr
#'                       ================= \tab === \tab ===========================================                                   \cr
#'          }
#'        }
#'        \item{}{
#'          \tabular{rcl}{
#'            \bold{\emph{DCCOverride}} =  \tab     \tab                                                             \cr
#'                                       1 \tab     \tab if \code{DCC} is missing or NA or not element of c(1:16).   \cr
#'                                         \tab     \tab ________________________________________________            \cr
#'                                         \tab     \tab \emph{Note:}                                                \cr
#'                                         \tab     \tab If the program cannot process the provided day count        \cr
#'                                         \tab     \tab identifier \code{DCC}, it overrides it with \code{DCC} = 2. \cr
#'                                         \tab     \tab ________________________________________________            \cr
#'                                         \tab     \tab \emph{Output:}                                              \cr
#'                                         \tab     \tab as if \code{DCC} = 2 was provided initially.                \cr
#'                     ___________________ \tab ___ \tab ________________________________________________            \cr
#'                                       0 \tab     \tab , else.                                                     \cr
#'                       ================= \tab === \tab ===========================================                 \cr
#'          }
#'        }
#'        \item{}{
#'          \tabular{rcl}{
#'            \bold{\emph{NoCoups}} =  \tab     \tab                                                                      \cr
#'                                   1 \tab     \tab , if there are no coupon payments between the provided               \cr
#'                                     \tab     \tab issue date (\code{Em}) and the maturity date (\code{Mat}), but the   \cr
#'                                     \tab     \tab provided (\code{CpY}) is not zero.                                   \cr
#'                                     \tab     \tab ________________________________________________                     \cr
#'                                     \tab     \tab \emph{Output:}                                                       \cr
#'                                     \tab     \tab \emph{RealDates} \code{= (Em,Mat)}, \emph{CoupDates} \code{= (Mat)}, \cr
#'                                     \tab     \tab \emph{AnnivDates} contains \code{Mat} and has either                 \cr
#'                                     \tab     \tab \code{2} or \code{3} elements, \emph{FCPType = LCPType} and          \cr
#'                                     \tab     \tab can be \code{"short"}, \code{"regular"} or \code{"long"}.            \cr
#'                 ___________________ \tab ___ \tab ________________________________________________                     \cr
#'                                   0 \tab     \tab , else.                                                              \cr
#'                   ================= \tab === \tab ===========================================                          \cr
#'         }
#'       }
#'     }
#'   }
#' }
#'
#' @references
#' \enumerate{
#'   \item{Djatschenko, Wadim, The Nitty Gritty of Bond Valuation: A Generalized Methodology for Fixed Coupon Bond Analysis Allowing for Irregular Periods and Various Day Count Conventions (November 5, 2018). Available at SSRN: https://ssrn.com/abstract=3205167.}
#' }
#'
#' @examples
#' data(SomeBonds2016)
#'
#' # Applying the function AnnivDates to the data frame SomeBonds2016.
#' system.time(
#'   FullAnalysis<-apply(SomeBonds2016[,c('Issue.Date','Mat.Date','CpY.Input','FIPD.Input',
#'   'LIPD.Input','FIAD.Input','RV.Input','Coup.Input','DCC.Input','EOM.Input')],1,function(y)
#'   AnnivDates(y[1],y[2],y[3],y[4],y[5],y[6],y[7],y[8],y[9],y[10],RegCF.equal=1)),
#' gcFirst = TRUE)
#' # warnings are due to apply's conversion of the variables' classes in
#' # SomeBonds2016 to class "character"
#'
#' # The output stored in FullAnalysis ist a nested list.
#' # Lets look at what is stored in FullAnalysis for a random bond:
#' randombond<-sample(c(1:nrow(SomeBonds2016)),1)
#' FullAnalysis[[randombond]]
#'
#' # Extracting the data frame Warnings:
#' AllWarnings<-do.call(rbind,lapply(FullAnalysis, `[[`, 1))
#' summary(AllWarnings)
#' # binding the Warnings to the bonds
#' BondsWithWarnings<-cbind(SomeBonds2016,AllWarnings)
#'
#' # Extracting the data frame Traits:
#' AllTraits<-do.call(rbind,lapply(FullAnalysis, `[[`, 2))
#' summary(AllTraits)
#' # binding the Traits to the bonds
#' BondsWithTraits<-cbind(SomeBonds2016,AllTraits)
#'
#' # Extracting the data frame AnnivDates:
#' AnnivDates<-lapply(lapply(FullAnalysis, `[[`, 3), `[[`, 5)
#' AnnivDates<-lapply(AnnivDates, `length<-`, max(lengths(AnnivDates)))
#' AnnivDates<-as.data.frame(do.call(rbind, AnnivDates))
#' AnnivDates<-as.data.frame(lapply(AnnivDates, as.Date, as.Date(AllTraits$DateOrigin[1])))
#' # binding the AnnivDates to the bonds:
#' BondsWithAnnivDates<-cbind(SomeBonds2016,AnnivDates)
#'
#' # Extracting the data frames PaySched for each bond and creating a panel:
#' CoupSched<-lapply(FullAnalysis, `[[`, 4)
#' CoupSchedPanel<-SomeBonds2016[rep(row.names(SomeBonds2016),sapply(CoupSched, nrow)),]
#' CoupSched<-as.data.frame(do.call(rbind, CoupSched))
#' CoupSchedPanel<-cbind(CoupSchedPanel,CoupSched)
#'
#'
#' @import Rcpp
#' @importFrom stats na.omit
#' @import timeDate
#' @export
#' @useDynLib BondValuation, .registration = TRUE
#'
AnnivDates<-function(Em=as.Date(NA),Mat=as.Date(NA),CpY=as.numeric(NA),FIPD=as.Date(NA),LIPD=as.Date(NA),FIAD=as.Date(NA),RV=as.numeric(NA),Coup=as.numeric(NA),DCC=as.numeric(NA),EOM=as.numeric(NA),DateOrigin=as.Date("1970-01-01"),InputCheck=1,FindEOM=FALSE,RegCF.equal=0) {
  if (length(Em)>1) {
    arglist<-Em
    argnames<-c("Em","Mat","CpY","FIPD","LIPD","FIAD","RV","Coup","DCC","EOM","DateOrigin","InputCheck","FindEOM","RegCF.equal")
    for (i in c(1:length(arglist))) {
      assign(argnames[i],arglist[i])
    }
  }
  # Checking whether the arguments are provided as the desired classes.
  # This is necessary to ensure that the function can be run with apply(), which transforms all input to
  # class "character" before it starts calcuations.
  # The code can deal with arguments provided as classes character, numeric or Date with format "%Y-%m-%d".
  # Otherwise the execution is cancelled with the appropriate error message.
  #
  if (InputCheck==1) {
    CheckedInput<-InputFormatCheck(Em=Em,Mat=Mat,CpY=CpY,FIPD=FIPD,LIPD=LIPD,FIAD=FIAD,RV=RV,Coup=Coup,DCC=DCC,EOM=EOM,DateOrigin=DateOrigin)
    Em<-CheckedInput$Em
    Mat<-CheckedInput$Mat
    CpY<-CheckedInput$CpY
    FIPD<-CheckedInput$FIPD
    LIPD<-CheckedInput$LIPD
    FIAD<-CheckedInput$FIAD
    RV<-CheckedInput$RV
    Coup<-CheckedInput$Coup
    DCC<-CheckedInput$DCC
    EOM<-CheckedInput$EOM
    DateOrigin<-CheckedInput$DateOrigin
  }
  Em_Orig<-Em
  FIPD_Orig<-FIPD
  LIPD_Orig<-LIPD
  EOM_Orig<-EOM
  est_EOM<-as.numeric(NA)
  EOM_used<-as.numeric(NA)
  Em_FIAD_differ<-0
  EmMatMissing<-0
  CpYOverride<-0
  RV_set100percent<-0
  NegLifeFlag<-0
  ZeroFlag<-0
  Em_Mat_SameMY<-0
  ChronErrorFlag<-0
  FIPD_LIPD_equal<-0   # if this is 1, it holds Em < FIPD = LIPD < Mat;
                       # Em = FIPD = LIPD and FIPD = LIPD = Mat are not considered here
  IPD_CpY_Corrupt<-0
  EOMOverride<-0
  EOM_Deviation<-0
  DCCOverride<-0
  NoCoups<-0
  est_FIPD<-as.Date(NA)
  est_LIPD<-as.Date(NA)
  Refer<-as.Date(NA)
  PaySched<-as.numeric(NA)
  AnnivDates<-as.Date(NA)
  CoupDates<-as.Date(NA)
  RealDates<-as.Date(NA)
  AD_indexes<-as.numeric(NA)
  RD_indexes<-as.numeric(NA)
  CD_indexes<-as.numeric(NA)
  index_FIPD<-as.numeric(NA)
  index_LIPD<-as.numeric(NA)
  FCPType<-as.character(NA)
  LCPType<-as.character(NA)
  f<-as.numeric(NA)
  l<-as.numeric(NA)
  if ((!(is.na(FIAD)))&(!(FIAD==Em))) {
    Em<-FIAD
    Em_FIAD_differ<-1
  }
  if ((missing(CpY))|(is.na(CpY))) {
    CpY<-2
    CpYOverride<-1
    warning("Number of interest payments per year (CpY) is missing or NA. CpY is set 2!")
  } else {
    if (!(is.element(CpY,c(0,1,2,3,4,6,12)))) {
      CpY<-2
      CpYOverride<-1
      warning("Number of interest payments per year (CpY) is not element of {0,1,2,3,4,6,12}!
              CpY is set 2!")
    }
  }
  if ((missing(RV))|(is.na(RV))) {
    RV<-100
    RV_set100percent<-1
    warning("Redemption value (RV) is missing or NA. RV is set 100!")
  }
  if (CpY==0) {
    Coup<-as.numeric(0)
    CpY<-as.numeric(1)
    # FIPD<-as.Date(NA)
    # LIPD<-as.Date(NA)
    # AnnivDates<-c(Em,Mat)
    # AnnivDates<-AnnivDates[!duplicated(AnnivDates)]
    # CoupDates<-c(Mat)
    # CoupDates<-CoupDates[!duplicated(CoupDates)]
    # RealDates<-c(Em,Mat)
    # RealDates<-RealDates[!duplicated(RealDates)]
    ZeroFlag<-1
    warning("This is a Zero Coupon bond!
            # CoupDates are (Mat), RealDates are (Em,Mat), AnnivDates are (Em,Mat).
            # The types of the first and last coupon periods are assigned NA!"
    )
  }
  if ((missing(Em))|(is.na(Em))) {
    EmMatMissing<-1
    warning("Issue date (Em) is missing or NA. NA created!")
  } else {
    if ((missing(Mat))|(is.na(Mat))) {
      EmMatMissing<-1
      warning("Maturity date (Mat) is missing or NA. NA created!")
    } else {
      if ((!(Mat>Em))) {
        NegLifeFlag<-1
        warning("Issue date (Em) is not before maturity date (Mat)! NA created!")
      } else {
        AtomVector_Em<-as.numeric(unlist(strsplit(as.character(Em),split = "-")))
        Atom1Em<-AtomVector_Em[1]
        Atom2Em<-AtomVector_Em[2]
        Atom3Em<-AtomVector_Em[3]
        EmHelp<-as.Date(paste(Atom1Em,Atom2Em,1,sep="-"))
        LDM_Em<-as.numeric(Date_LDM(c(Atom1Em,Atom2Em,Atom3Em)))
        LDM_Em<-as.Date(paste(LDM_Em[1],LDM_Em[2],LDM_Em[3],sep="-"))
        AtomVector_Mat<-as.numeric(unlist(strsplit(as.character(Mat),split = "-")))
        Atom1Mat<-AtomVector_Mat[1]
        Atom2Mat<-AtomVector_Mat[2]
        Atom3Mat<-AtomVector_Mat[3]
        MatHelp<-as.Date(paste(Atom1Mat,Atom2Mat,15,sep="-"))
        LDM_Mat<-as.numeric(Date_LDM(c(Atom1Mat,Atom2Mat,Atom3Mat)))
        LDM_Mat<-as.Date(paste(LDM_Mat[1],LDM_Mat[2],LDM_Mat[3],sep="-"))
        if ((Atom1Em==Atom1Mat)&(Atom2Em==Atom2Mat)) {
          Em_Mat_SameMY<-1
          warning("The issue date and the maturity date are in the same month of the same year!
                  CoupDates are (Mat), RealDates are (Em,Mat). The types of the first and last
                  coupon periods are assigned \"short\"!")
        }
        # N_months is essential for further calculations
        if (CpY==1) {N_months<-12}
        if (CpY==2) {N_months<-6}
        if (CpY==3) {N_months<-4}
        if (CpY==4) {N_months<-3}
        if (CpY==6) {N_months<-2}
        if (CpY==12) {N_months<-1}
        # the following code tests for the correct chronological order of the provided dates
        ChronoVec_input<-c(Em,FIPD,LIPD,Mat)
        ChronoVec_input<-na.omit(ChronoVec_input)
        ChronoVec_test<-sort(ChronoVec_input)
        if (any(ChronoVec_input!=ChronoVec_test)) {
          FIPD<-as.Date(NA)
          LIPD<-as.Date(NA)
          ChronErrorFlag<-1
          warning("The date inputs are in a wrong chronological order!
                  FIPD and LIPD dropped.
                  Note: The correct ascending chronological order is
                  issue date (Em), first interest payment date (FIPD), last interest payment date (LIPD), maturity date (Mat).
                  (Please note, that issue date (Em) is substituted by the first interest accrual date (FIAD)
                  if FIAD is available and different from Em.)")
        }
        # if FIPD and/or LIPD are not dropped, it now holds Em <= FIPD <= LIPD <= Mat and Em < Mat
        # the following code deals with special cases of the location of the dates
        if (is.na(FIPD)) {
          if (!(is.na(LIPD))) {     # Case 1: LIPD is available but FIPD is not
            if (!(Em<LIPD)) {
              LIPD<-as.Date(NA)      # Em is not prior to LIPD --> LIPD is assigned NA
            } else {
              if (!(LIPD<Mat)) {
                LIPD<-as.Date(NA)    # LIPD is not prior to Mat --> LIPD is assigned NA
              }
            }
          }
        }
        if (!(is.na(FIPD))) {
          if (is.na(LIPD)) {        # Case 2: FIPD is available but LIPD is not
            if (!(Em<FIPD)) {
              FIPD<-as.Date(NA)      # Em is not prior to FIPD --> FIPD is assigned NA
            } else {
              if (!(FIPD<Mat)) {
                FIPD<-as.Date(NA)    # FIPD is not prior to Mat --> FIPD is assigned NA
              }
            }
          }
        }
        if ((!is.na(FIPD))&(!is.na(LIPD))) {     # Case 3: FIPD and LIPD are both available
          if (Em<FIPD) {
            if (FIPD<LIPD) {
              if (LIPD<Mat) {
                # it holds Em < FIPD < LIPD < Mat --> do nothing
              } else { # i.e. if (LIPD==Mat)
                # it holds Em < FIPD < LIPD = Mat --> LIPD is redundant and is dropped
                LIPD<-as.Date(NA)
              }
            } else { # i.e. if (FIPD==LIPD)
              if (LIPD<Mat) {
                # it holds Em < FIPD = LIPD < Mat --> this case is considered further below,
                # here only the flag is changed
                FIPD_LIPD_equal<-1
              } else { # i.e. if (LIPD==Mat)
                # it holds Em < FIPD = LIPD = Mat --> FIPD and LIPD are redundant and dropped
                FIPD<-as.Date(NA)
                LIPD<-as.Date(NA)
              }
            }
          } else { # i.e. if (Em==FIPD)
            if (FIPD<LIPD) {
              if (LIPD<Mat) {
                # it holds Em = FIPD < LIPD < Mat --> LIPD is redundant and is dropped
                FIPD<-as.Date(NA)
              } else { # i.e. if (LIPD==Mat)
                # it holds Em = FIPD < LIPD = Mat --> FIPD and LIPD are redundant and dropped
                FIPD<-as.Date(NA)
                LIPD<-as.Date(NA)
              }
            } else { # i.e. if (FIPD==LIPD)
              # it holds Em = FIPD = LIPD < Mat --> FIPD and LIPD are redundant and dropped
              FIPD<-as.Date(NA)
              LIPD<-as.Date(NA)
            }
          }
        }
        # If FIPD and LIPD are both dropped, it now holds:
        # Em < Mat and FIPD_LIPD_equal = 0
        # If either FIPD or LIPD is dropped, it now holds:
        # Em < FIPD (or LIPD) < Mat and FIPD_LIPD_equal = 0
        # If FIPD and LIPD are not both dropped, two cases can occur:
        # 1. Em < FIPD = LIPD < Mat and FIPD_LIPD_equal = 1
        # 2. Em < FIPD < LIPD < Mat and FIPD_LIPD_equal = 0
        #
        # But still it can happen, that FIPD and LIPD are inconsistent with CpY, such that the assumption,
        # that the anniversary dates between FIPD and LIPD are evenly distributed, is violated.
        # This happens,
        # 1. if FIPD and LIPD are in the same month of the same year but not on the same day, or
        # 2. if the month difference between FIPD and LIPD is not a multiple of the number of months implied by CpY (N_months), or
        # 3. if FIPD and LIPD are not both last day in month but their %d differ and the day figure difference between FIPD and LIPD is not due to different month lengths.
        #
        if (!is.na(FIPD)) {
          AtomVector_FIPD<-as.numeric(unlist(strsplit(as.character(FIPD),split = "-")))
          Atom1FIPD<-AtomVector_FIPD[1]
          Atom2FIPD<-AtomVector_FIPD[2]
          Atom3FIPD<-AtomVector_FIPD[3]
          FIPDHelp<-as.Date(paste(Atom1FIPD,Atom2FIPD,15,sep="-"))
          LDM_FIPD<-as.numeric(Date_LDM(c(Atom1FIPD,Atom2FIPD,Atom3FIPD)))
          LDM_FIPD<-as.Date(paste(LDM_FIPD[1],LDM_FIPD[2],LDM_FIPD[3],sep="-"))
        }
        if (!is.na(LIPD)) {
          AtomVector_LIPD<-as.numeric(unlist(strsplit(as.character(LIPD),split = "-")))
          Atom1LIPD<-AtomVector_LIPD[1]
          Atom2LIPD<-AtomVector_LIPD[2]
          Atom3LIPD<-AtomVector_LIPD[3]
          LIPDHelp<-as.Date(paste(Atom1LIPD,Atom2LIPD,15,sep="-"))
          LDM_LIPD<-as.numeric(Date_LDM(c(Atom1LIPD,Atom2LIPD,Atom3LIPD)))
          LDM_LIPD<-as.Date(paste(LDM_LIPD[1],LDM_LIPD[2],LDM_LIPD[3],sep="-"))
        }
        if ((!is.na(FIPD))&(!is.na(LIPD))) {
          if ((Atom1FIPD==Atom1LIPD)&(Atom2FIPD==Atom2LIPD)&(Atom3FIPD!=Atom3LIPD)) {
            FIPD<-as.Date(NA)
            LIPD<-as.Date(NA)
            IPD_CpY_Corrupt<-1
            warning("FIPD and LIPD are inconsistent with CpY! FIPD and LIPD are dropped.
                    Note: The assumption, that the anniversary dates between FIPD and LIPD are evenly distributed, is violated!
                    Cause: FIPD and LIPD are in the same month of the same year but not on the same day.")
          } else {
            if (abs((as.POSIXlt(LIPD)$year*12+as.POSIXlt(LIPD)$mon+1)-(as.POSIXlt(FIPD)$year*12+as.POSIXlt(FIPD)$mon+1))%%N_months!=0) {
              FIPD<-as.Date(NA)
              LIPD<-as.Date(NA)
              IPD_CpY_Corrupt<-1
              warning("FIPD and LIPD are inconsistent with CpY! FIPD and LIPD are dropped.
                      Note: The assumption, that the anniversary dates between FIPD and LIPD are evenly distributed, is violated!
                      Cause: The month difference between FIPD and LIPD is not a multiple of the number of months implied by CpY.")
            } else {
              if (Atom3FIPD==Atom3LIPD) {   # i.e. FIPD and LIPD have the same %d
                # do nothing
              } else {   # i.e. FIPD and LIPD have different %d
                if (LDM_FIPD==FIPD) {   # i.e. FIPD is last day in month
                  if (LDM_LIPD==LIPD) {   # i.e. LIPD is last day in month
                    # do nothing
                  } else {   # i.e. LIPD is not last day in month
                    if ((Atom3FIPD+1)==Atom3LIPD) {   # i.e. Increasing the %d of FIPD by 1 results in the %d of LIPD.
                      # do nothing
                    } else {   # i.e. Increasing the %d of FIPD by 1 does not result in the %d of LIPD.
                      if ((Atom3FIPD+2)==Atom3LIPD) {   # i.e. Increasing the %d of FIPD by 2 results in the %d of LIPD.
                        # do nothing
                      } else {   # i.e. Increasing the %d of FIPD by 2 does not result in the %d of LIPD.
                        FIPD<-as.Date(NA)
                        LIPD<-as.Date(NA)
                        IPD_CpY_Corrupt<-1
                        warning("FIPD and LIPD are inconsistent with CpY! FIPD and LIPD are dropped.
                                Note: The assumption, that the anniversary dates between FIPD and LIPD are evenly distributed, is violated!
                                Cause: FIPD and LIPD have different %d. FIPD is last day in month. LIPD is not last day in month.
                                The day figure difference between FIPD and LIPD is not due to different month lengths.")
                      }
                      }
                      }
                      } else {   # i.e. FIPD is not last day in month
                        if (LDM_LIPD==LIPD) {   # i.e. LIPD is last day in month
                          if ((Atom3LIPD+1)==Atom3FIPD) {   # i.e. Increasing the %d of LIPD by 1 results in the %d of FIPD.
                            # do nothing
                          } else {   # i.e. Increasing the %d of LIPD by 1 does not result in the %d of FIPD.
                            if ((Atom3LIPD+2)==Atom3FIPD) {   # i.e. Increasing the %d of LIPD by 2 results in the %d of FIPD.
                              # do nothing
                            } else {   # i.e. Increasing the %d of LIPD by 2 does not result in the %d of FIPD.
                              FIPD<-as.Date(NA)
                              LIPD<-as.Date(NA)
                              IPD_CpY_Corrupt<-1
                              warning("FIPD and LIPD are inconsistent with CpY! FIPD and LIPD are dropped.
                                      Note: The assumption, that the anniversary dates between FIPD and LIPD are evenly distributed, is violated!
                                      Cause: FIPD and LIPD have different %d. FIPD is not last day in month. LIPD is last day in month.
                                      The day figure difference between FIPD and LIPD is not due to different month lengths.")
                            }
                            }
                            } else {   # i.e. LIPD is not last day in month
                              FIPD<-as.Date(NA)
                              LIPD<-as.Date(NA)
                              IPD_CpY_Corrupt<-1
                              warning("FIPD and LIPD are inconsistent with CpY! FIPD and LIPD are dropped.
                                      Note: The assumption, that the anniversary dates between FIPD and LIPD are evenly distributed, is violated!
                                      Cause: FIPD and LIPD have different %d. FIPD is not last day in month. LIPD is not last day in month.")
                            }
                          }
                        }
                    }
              }
        }
        # The following code determines est_EOM based on the available calendar dates.
        if (is.na(FIPD)) {
          if (is.na(LIPD)) {      # Case 1: FIPD and LIPD are both NA
            if (LDM_Mat==Mat) {
              est_EOM<-1
            } else {
              est_EOM<-0
            }
          } else {                # Case 2: FIPD is NA and LIPD is available
            if (LDM_LIPD==LIPD) {
              est_EOM<-1
            } else {
              est_EOM<-0
            }
          }
        } else {
          if (is.na(LIPD)) {      # Case 3: FIPD is available and LIPD is NA
            if (LDM_FIPD==FIPD) {
              est_EOM<-1
            } else {
              est_EOM<-0
            }
          } else {                # Case 4: FIPD and LIPD are both available
            if ((LDM_FIPD==FIPD)&(LDM_LIPD==LIPD)) {
              est_EOM<-1
            } else {
              est_EOM<-0
            }
          }
        }
        if (is.na(EOM_Orig)) {
          EOMOverride<-1
          EOM<-est_EOM
          warning(paste("EOM was not provided or NA! EOM is set",EOM,".
                        Note: The available calandar dates suggest that EOM =",est_EOM,"."))
        } else {
          if (EOM_Orig!=est_EOM) {
            EOM_Deviation<-1
            if (FindEOM==TRUE) {
              EOM<-est_EOM
              EOMOverride<-1
              warning(paste("The available calandar dates suggest that EOM =",est_EOM,".
                            Option FindEOM = TRUE is active. EOM is set",est_EOM,"."))
            } else {
              if ((est_EOM==0)&(EOM==1)) {
                EOM<-est_EOM
                EOMOverride<-1
                warning(paste("The provided EOM =",EOM_Orig,"conflicts with one or several provided calendar dates.
                              EOM is set",est_EOM,"."))
              }
              if ((est_EOM==1)&(EOM==0)) {
                warning(paste("The available calandar dates suggest that EOM =",est_EOM,".
                              Option FindEOM = FALSE is active. Provided EOM is not overridden and remains EOM =",EOM,"."))
              }
              }
              }
        }
        # If DCC is not provided or NA or not element of {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16}, the following code sets it 2 (Act/Act (ICMA)).
        if ((missing(DCC))|(is.na(DCC))) {
          DCC<-2
          DCCOverride<-1
          warning("The day count indentifier (DCC) is missing or NA. DCC is set 2 (Act/Act (ICMA))!")
        } else {
          if (!(is.element(DCC,c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)))) {
            DCC<-2
            DCCOverride<-1
            warning("The day count indentifier (DCC) is not element of {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16}!
                    DCC is set 2 (Act/Act (ICMA))!")
          }
          }
        ######>>>>>
        ######>>>>> here the computation of the date vectors begins
        ######>>>>>
        if (is.na(FIPD)) { # i.e. if FIPD is not given
          if (is.na(LIPD)) { # i.e. if LIPD is not given
            # If neither FIPD nor LIPD are available the following code evaluates the bond based only upon the "essential" variables Em, Mat, CpY.
            # It creates:
            # estimated FIPD and estimated LIPD
            # Vector RealDates: (Em,estimated FIPD,regular coupon payment dates,estimated LIPD,Mat)
            # Vector CoupDates: (estimated FIPD,regular coupon payment dates,estimated LIPD,Mat)
            # Vector AnnivDates: (notional Em, coupon payment dates incl. estimated FIPD and estimated LIPD, Mat)
            #                    Em is only element of AnnivDates, if the bond has a regular first coupon period.
            # numeric E: the issue date's position in time relative to the bond's notional issue date
            # numeric M: the maturity date's position in time relative to the bond's notional issue date
            # numeric f: length of the first coupon period as a fraction of a regular coupon period
            # numeric l: length of the final coupon period as a fraction of a regular coupon period
            #
            # Assumptions:
            # >>> All coupon periods but the first are "regular". <<<
            # >>> The value of EOM determines the location of all coupon payment dates. <<<
            # >>> The interest payment dates between the estimated FIPD and Mat are evenly distributed. <<<
            # >>> Since FIPD is not given, it is impossible to distinguish between a "short" and "long" odd first coupon period,
            # without an assumption on the number of interest payment dates.
            # Consequently the first coupon period is either "regular" or "short". <<<
            #
            ###
            ### Creating AnnivDates
            ###
            AnnivDates<-rev(seq(MatHelp,EmHelp,by=paste("-",N_months," months",sep="")))
            AnnivDates<-sort(na.omit(AnnivDates[!duplicated(AnnivDates)]))
            # assigning the reference date that determines the day figures of all AnnivDates
            Atom1Refer<-Atom1Mat
            Atom2Refer<-Atom2Mat
            Atom3Refer<-Atom3Mat
            if (EOM==1) {
              AnnivDates<-as.Date(timeLastDayInMonth(AnnivDates))
            } else {
              AnnivDates_A<-as.Date(ISOdatetime(t(atoms(as.timeDate(AnnivDates))[1]),t(atoms(as.timeDate(AnnivDates))[2]),Atom3Refer,12,0,0))
              if (length(which(is.na(AnnivDates_A)))!=0) {
                # Assuming %d of Mat as the %d of all AnnivDates produced NAs. NAs are substituted by the respective last day in month.
                nas<-which(is.na(AnnivDates_A))
                AnnivDates_B<-as.Date(timeLastDayInMonth(AnnivDates[nas]))
                AnnivDates<-sort(na.omit(append(AnnivDates_A,AnnivDates_B)))
              } else {
                AnnivDates<-AnnivDates_A
              }
            }
            AnnivDates<-sort(na.omit(AnnivDates[!duplicated(AnnivDates)]))
            # creating the anniversary date preceding Em
            if (AnnivDates[1]!=Em) {
              AtomVector_AD1<-as.numeric(unlist(strsplit(as.character(AnnivDates[1]),split = "-")))
              Atom1AD1<-AtomVector_AD1[1]
              Atom2AD1<-AtomVector_AD1[2]
              Atom3AD1<-AtomVector_AD1[3]
              PrevDate<-as.numeric(CppPrevDate(c(Atom1AD1,Atom2AD1,Atom3AD1,Atom1Em,Atom2Em,Atom3Em,Atom1Refer,Atom2Refer,Atom3Refer,CpY,EOM)))
              if (any(is.na(PrevDate))) {
                PrevDate<-as.Date(NA)
              } else {
                PrevDate<-as.Date(paste(PrevDate[1],PrevDate[2],PrevDate[3],sep="-"))
              }
              AnnivDates<-c(PrevDate,AnnivDates)
              AnnivDates<-sort(na.omit(AnnivDates[!duplicated(AnnivDates)]))   # final AnnivDates created!
            }
            AD_indexes<-c(1:length(AnnivDates))-1
            AD_List<-list(AnnivDates,AD_indexes)
            index_FIPD<-1
            index_LIPD<-AD_indexes[length(AD_indexes)-1]
            ###
            ### now estimating FIPD and LIPD and creating the vectors RealDates and CoupDates
            ###
            if (length(AnnivDates)==2) {  # this is the case when there are no interest payments between Em and Mat
              NoCoups<-1
              est_FIPD<-as.Date(NA)
              est_LIPD<-as.Date(NA)
              warning(paste("There are no interest payments between the issue date ( Em =", Em,") and the maturity date ( Mat =",Mat,"), but the number of interest payments per year is CpY =", CpY,". The estimated FIPD and LIPD are assigned NA."))
            } else {
              est_FIPD<-AnnivDates[2]
              est_LIPD<-AnnivDates[length(AnnivDates)-1]
            }
            AD_in_RD<-as.numeric(which((!(AnnivDates<est_FIPD))&(!(AnnivDates>est_LIPD))))
            if (length(length(AD_in_RD))==0) {
              RealDates<-c(Em,Mat)
            } else {
              RealDates<-c(Em,AnnivDates[AD_in_RD],Mat)
            }
            RealDates<-RealDates[!duplicated(RealDates)]
            CoupDates<-RealDates[-1]
          } else { # i.e. if LIPD is given
            # If LIPD is available but FIPD is not, the following code creates:
            # estimated FIPD
            # Vector RealDates: (Em,estimated FIPD,regular coupon payment dates,LIPD,Mat)
            # Vector CoupDates: (estimated FIPD,regular coupon payment dates,LIPD,Mat)
            # Vector AnnivDates: (notional Em, coupon payment dates incl. estimated FIPD and LIPD,
            #                    possibly notional coupon payment dates between LIPD and Mat,
            #                    notional Mat)
            #                    Em is only element of AnnivDates, if the bond has a regular first coupon period.
            #                      If Em is element of AnnivDates, it is the first element.
            #                    Mat is element of AnnivDates, if the bond has a regular final coupon period or
            #                      if the bond has a long final coupon period whose length is a multiple of a
            #                      regular coupon period.
            #                      If Mat is element of AnnivDates, it is the last element.
            # numeric E: the issue date's position in time relative to the bond's notional issue date
            # numeric M: the maturity date's position in time relative to the bond's notional issue date
            # numeric f: length of the first coupon period as a fraction of a regular coupon period
            # numeric l: length of the final coupon period as a fraction of a regular coupon period
            #
            # Since LIPD is available, the length of the last coupon period is determined by LIPD and Mat,
            # so it can be "long", "regular" or "short".
            #
            # To find the interest payment dates between the estimated FIPD and LIPD the following assumptions are met:
            # >>> The interest payment dates between the estimated FIPD and LIPD are evenly distributed. <<<
            # >>> The value of EOM determines the location of all coupon payment dates. <<<
            # >>> Since FIPD is not given, it is impossible to distinguish between a "short" and "long" odd first coupon period,
            # without an assumption on the number of interest payment dates.
            # Consequently the first coupon period is either "regular" or "short". <<<
            #
            ###
            ### Creating AnnivDates
            ###
            AnnivDates_Head<-rev(seq(LIPDHelp,EmHelp,by=paste("-",N_months," months",sep="")))
            AnnivDates_Tail<-seq(LIPDHelp,LDM_Mat,by=paste(N_months," months",sep=""))
            AnnivDates<-c(AnnivDates_Head,LIPD,AnnivDates_Tail)
            AnnivDates<-sort(na.omit(AnnivDates[!duplicated(AnnivDates)]))
            # assigning the reference date that determines the day figures of all AnnivDates
            Atom1Refer<-Atom1LIPD
            Atom2Refer<-Atom2LIPD
            Atom3Refer<-Atom3LIPD
            if (EOM==1) {
              AnnivDates<-as.Date(timeLastDayInMonth(AnnivDates))
            } else {
              AnnivDates_A<-as.Date(ISOdatetime(t(atoms(as.timeDate(AnnivDates))[1]),t(atoms(as.timeDate(AnnivDates))[2]),Atom3Refer,12,0,0))
              if (length(which(is.na(AnnivDates_A)))!=0) {
                # Assuming %d of LIPD as the %d of all AnnivDates produced NAs. NAs are substituted by the respective last day in month.
                nas<-which(is.na(AnnivDates_A))
                AnnivDates_B<-as.Date(timeLastDayInMonth(AnnivDates[nas]))
                AnnivDates<-sort(na.omit(append(AnnivDates_A,AnnivDates_B)))
              } else {
                AnnivDates<-AnnivDates_A
              }
            }
            AnnivDates<-sort(AnnivDates[!duplicated(AnnivDates)])
            # creating the anniversary date preceding Em
            if (AnnivDates[1]!=Em) {
              AtomVector_AD1<-as.numeric(unlist(strsplit(as.character(AnnivDates[1]),split = "-")))
              Atom1AD1<-AtomVector_AD1[1]
              Atom2AD1<-AtomVector_AD1[2]
              Atom3AD1<-AtomVector_AD1[3]
              PrevDate<-as.numeric(CppPrevDate(c(Atom1AD1,Atom2AD1,Atom3AD1,Atom1Em,Atom2Em,Atom3Em,Atom1Refer,Atom2Refer,Atom3Refer,CpY,EOM)))
              if (any(is.na(PrevDate))) {
                PrevDate<-as.Date(NA)
              } else {
                PrevDate<-as.Date(paste(PrevDate[1],PrevDate[2],PrevDate[3],sep="-"))
              }
              AnnivDates<-c(PrevDate,AnnivDates)
              AnnivDates<-sort(na.omit(AnnivDates[!duplicated(AnnivDates)]))
            }
            # creating the anniversary date succeeding Mat
            if (AnnivDates[length(AnnivDates)]!=Mat) {
              AtomVector_ADfin<-as.numeric(unlist(strsplit(as.character(AnnivDates[length(AnnivDates)]),split = "-")))
              Atom1ADfin<-AtomVector_ADfin[1]
              Atom2ADfin<-AtomVector_ADfin[2]
              Atom3ADfin<-AtomVector_ADfin[3]
              SuccDate<-as.numeric(CppSuccDate(c(Atom1ADfin,Atom2ADfin,Atom3ADfin,Atom1Mat,Atom2Mat,Atom3Mat,Atom1Refer,Atom2Refer,Atom3Refer,CpY,EOM)))
              if (any(is.na(SuccDate))) {
                SuccDate<-as.Date(NA)
              } else {
                SuccDate<-as.Date(paste(SuccDate[1],SuccDate[2],SuccDate[3],sep="-"))
              }
              AnnivDates<-c(AnnivDates,SuccDate)
              AnnivDates<-sort(na.omit(AnnivDates[!duplicated(AnnivDates)]))   # final AnnivDates created!
            }
            AD_indexes<-c(1:length(AnnivDates))-1
            AD_List<-list(AnnivDates,AD_indexes)
            index_FIPD<-1
            index_LIPD<-AD_List[[2]][which(AD_List[[1]]==LIPD)]
            ###
            ### now estimating FIPD and creating the vectors RealDates and CoupDates
            ###
            #
            # AnnivDates has at least 3 elements
            # AD1 can be Em, AD_Fin can be Mat, AD_PreFin can be LIPD
            #
            est_FIPD<-AnnivDates[2]
            AD_in_RD<-as.numeric(which((!(AnnivDates<est_FIPD))&(!(AnnivDates>LIPD))))
            if (length(length(AD_in_RD))==0) {
              RealDates<-c(Em,Mat)
            } else {
              RealDates<-c(Em,AnnivDates[AD_in_RD],Mat)
            }
            RealDates<-RealDates[!duplicated(RealDates)]
            CoupDates<-RealDates[-1]
          }
        } else { # i.e. if FIPD is given
          if (is.na(LIPD)) { # i.e. if LIPD is not given
            # If FIPD is available but LIPD is not, the following code creates:
            # estimated LIPD
            # Vector RealDates: (Em,FIPD,regular coupon payment dates,estimated LIPD,Mat)
            # Vector CoupDates: (FIPD,regular coupon payment dates,estimated LIPD,Mat)
            # Vector AnnivDates: (notional Em, possibly notional coupon payment dates between Em and FIPD,
            #                    coupon payment dates incl. FIPD and estimated LIPD, notional Mat)
            #                    Em is element of AnnivDates, if the bond has a regular first coupon period or
            #                      if the bond has a long first coupon period whose length is a multiple of a
            #                      regular coupon period.
            #                      If Em is element of AnnivDates, it is the first element.
            #                    Mat is only element of AnnivDates, if the bond has a regular final coupon period.
            #                      If Mat is element of AnnivDates, it is the last element.
            # numeric E: the issue date's position in time relative to the bond's notional issue date
            # numeric M: the maturity date's position in time relative to the bond's notional issue date
            # numeric f: length of the first coupon period as a fraction of a regular coupon period
            # numeric l: length of the final coupon period as a fraction of a regular coupon period
            #
            # Since FIPD is available, the length of the first coupon period is determined by Em and FIPD,
            # so it can be "long", "regular" or "short".
            #
            # To find the interest payment dates between FIPD and the estimated LIPD the following assumptions are met:
            # >>> The interest payment dates between FIPD and the estimated LIPD are evenly distributed. <<<
            # >>> The value of EOM determines the location of all coupon payment dates. <<<
            # >>> Since LIPD is not given, it is impossible to distinguish between a "short" and "long" odd final coupon period,
            # without an assumption on the number of interest payment dates.
            # Consequently the final coupon period is either "regular" or "short". <<<
            #
            ###
            ### Creating AnnivDates
            ###
            AnnivDates_Head<-rev(seq(FIPDHelp,EmHelp,by=paste("-",N_months," months",sep="")))
            AnnivDates_Tail<-seq(FIPDHelp,LDM_Mat,by=paste(N_months," months",sep=""))
            AnnivDates<-c(AnnivDates_Head,FIPD,AnnivDates_Tail)
            AnnivDates<-sort(na.omit(AnnivDates[!duplicated(AnnivDates)]))
            # assigning the reference date that determines the day figures of all AnnivDates
            Atom1Refer<-Atom1FIPD
            Atom2Refer<-Atom2FIPD
            Atom3Refer<-Atom3FIPD
            if (EOM==1) {
              AnnivDates<-as.Date(timeLastDayInMonth(AnnivDates))
            } else {
              AnnivDates_A<-as.Date(ISOdatetime(t(atoms(as.timeDate(AnnivDates))[1]),t(atoms(as.timeDate(AnnivDates))[2]),Atom3Refer,12,0,0))
              if (length(which(is.na(AnnivDates_A)))!=0) {
                # Assuming %d of LIPD as the %d of all AnnivDates produced NAs. NAs are substituted by the respective last day in month.
                nas<-which(is.na(AnnivDates_A))
                AnnivDates_B<-as.Date(timeLastDayInMonth(AnnivDates[nas]))
                AnnivDates<-sort(na.omit(append(AnnivDates_A,AnnivDates_B)))
              } else {
                AnnivDates<-AnnivDates_A
              }
            }
            AnnivDates<-sort(AnnivDates[!duplicated(AnnivDates)])
            # creating the anniversary date preceding Em
            if (AnnivDates[1]!=Em) {
              AtomVector_AD1<-as.numeric(unlist(strsplit(as.character(AnnivDates[1]),split = "-")))
              Atom1AD1<-AtomVector_AD1[1]
              Atom2AD1<-AtomVector_AD1[2]
              Atom3AD1<-AtomVector_AD1[3]
              PrevDate<-as.numeric(CppPrevDate(c(Atom1AD1,Atom2AD1,Atom3AD1,Atom1Em,Atom2Em,Atom3Em,Atom1Refer,Atom2Refer,Atom3Refer,CpY,EOM)))
              if (any(is.na(PrevDate))) {
                PrevDate<-as.Date(NA)
              } else {
                PrevDate<-as.Date(paste(PrevDate[1],PrevDate[2],PrevDate[3],sep="-"))
              }
              AnnivDates<-c(PrevDate,AnnivDates)
              AnnivDates<-sort(na.omit(AnnivDates[!duplicated(AnnivDates)]))
            }
            # creating the anniversary date succeeding Mat
            if (AnnivDates[length(AnnivDates)]!=Mat) {
              AtomVector_ADfin<-as.numeric(unlist(strsplit(as.character(AnnivDates[length(AnnivDates)]),split = "-")))
              Atom1ADfin<-AtomVector_ADfin[1]
              Atom2ADfin<-AtomVector_ADfin[2]
              Atom3ADfin<-AtomVector_ADfin[3]
              SuccDate<-as.numeric(CppSuccDate(c(Atom1ADfin,Atom2ADfin,Atom3ADfin,Atom1Mat,Atom2Mat,Atom3Mat,Atom1Refer,Atom2Refer,Atom3Refer,CpY,EOM)))
              if (any(is.na(SuccDate))) {
                SuccDate<-as.Date(NA)
              } else {
                SuccDate<-as.Date(paste(SuccDate[1],SuccDate[2],SuccDate[3],sep="-"))
              }
              AnnivDates<-c(AnnivDates,SuccDate)
              AnnivDates<-sort(na.omit(AnnivDates[!duplicated(AnnivDates)]))   # final AnnivDates created!
            }
            AD_indexes<-c(1:length(AnnivDates))-(which(AnnivDates==FIPD)-1)
            AD_List<-list(AnnivDates,AD_indexes)
            index_FIPD<-AD_List[[2]][which(AD_List[[1]]==FIPD)]
            index_LIPD<-AD_indexes[length(AD_indexes)-1]
            ###
            ### now estimating LIPD and creating the vectors RealDates and CoupDates
            ###
            #
            # AnnivDates has at least 3 elements
            # AD1 can be Em, AD_Fin can be Mat, AD_2 can be FIPD
            #
            est_LIPD<-AnnivDates[length(AnnivDates)-1]
            AD_in_RD<-as.numeric(which((!(AnnivDates<FIPD))&(!(AnnivDates>est_LIPD))))
            if (length(length(AD_in_RD))==0) {
              RealDates<-c(Em,Mat)
            } else {
              RealDates<-c(Em,AnnivDates[AD_in_RD],Mat)
            }
            RealDates<-RealDates[!duplicated(RealDates)]
            CoupDates<-RealDates[-1]
          } else { # i.e. if LIPD is given
            # If FIPD and LIPD are both available the following code creates:
            # Vector RealDates: (Em,FIPD,regular coupon payment dates,LIPD,Mat)
            # Vector CoupDates: (FIPD,regular coupon payment dates,LIPD,Mat)
            # Vector AnnivDates: (notional Em, possibly notional coupon payment dates between Em and FIPD,
            #                    coupon payment dates incl. FIPD and LIPD, possibly notional coupon payment
            #                    dates between LIPD and Mat, notional Mat)
            #                    Em is element of AnnivDates, if the bond has a regular first coupon period or
            #                      if the bond has a long first coupon period whose length is a multiple of a
            #                      regular coupon period.
            #                      If Em is element of AnnivDates, it is the first element.
            #                    Mat is element of AnnivDates, if the bond has a regular final coupon period or
            #                      if the bond has a long final coupon period whose length is a multiple of a
            #                      regular coupon period.
            #                      If Mat is element of AnnivDates, it is the last element.
            # numeric E: the issue date's position in time relative to the bond's notional issue date
            # numeric M: the maturity date's position in time relative to the bond's notional issue date
            # numeric f: length of the first coupon period as a fraction of a regular coupon period
            # numeric l: length of the final coupon period as a fraction of a regular coupon period
            #
            # Since FIPD and LIPD are both available, the lengths of the first and last coupon periods are determinate
            # and can be "long", "regular" or "short".
            #
            # To find the interest payment dates between FIPD and LIPD the following assumptions are met:
            # >>> The interest payment dates between FIPD and LIPD are evenly distributed. <<<
            # >>> The value of EOM determines the location of all coupon payment dates. <<<
            #
            # Notes:
            # If the %d of both FIPD and LIPD are equal or both last day in month there is no ambiguity on the concrete dates.
            # If the %d of FIPD and LIPD differ and are not both last day in month, and
            # the day figure difference between FIPD and LIPD is not due to different month lengths,
            # then the first assumption is violated and FIPD and LIPD are already dropped with IPD_CpY_Corrupt = 1.
            # If IPD_CpY_Corrupt = 0, then EOM = 0 and the higer day figure among FIPD and LIPD determines the
            # day figures of all anniversary dates.
            #
            #
            ###
            ### Creating AnnivDates
            ###
            AnnivDates_Head<-rev(seq(FIPDHelp,EmHelp,by=paste("-",N_months," months",sep="")))
            AnnivDates_Tail<-seq(FIPDHelp,LDM_Mat,by=paste(N_months," months",sep=""))
            AnnivDates<-c(AnnivDates_Head,FIPD,LIPD,AnnivDates_Tail)
            AnnivDates<-sort(na.omit(AnnivDates[!duplicated(AnnivDates)]))
            # assigning the reference date that determines the day figures of all AnnivDates
            if (EOM==0) {
              if ((LDM_FIPD==FIPD)&(LDM_LIPD==LIPD)) {
                if (Atom3LIPD>Atom3FIPD) {
                  Atom1Refer<-Atom1LIPD
                  Atom2Refer<-Atom2LIPD
                  Atom3Refer<-Atom3LIPD
                } else {
                  Atom1Refer<-Atom1FIPD
                  Atom2Refer<-Atom2FIPD
                  Atom3Refer<-Atom3FIPD
                }
              } else {
                if (LDM_LIPD==LIPD) {
                  Atom1Refer<-Atom1FIPD
                  Atom2Refer<-Atom2FIPD
                  Atom3Refer<-Atom3FIPD
                } else {
                  Atom1Refer<-Atom1LIPD
                  Atom2Refer<-Atom2LIPD
                  Atom3Refer<-Atom3LIPD
                }
              }
            } else {
              Atom1Refer<-Atom1LIPD
              Atom2Refer<-Atom2LIPD
              Atom3Refer<-Atom3LIPD
            }
            if (EOM==1) {
              AnnivDates<-as.Date(timeLastDayInMonth(AnnivDates))
            } else {
              AnnivDates_A<-as.Date(ISOdatetime(t(atoms(as.timeDate(AnnivDates))[1]),t(atoms(as.timeDate(AnnivDates))[2]),Atom3Refer,12,0,0))
              if (length(which(is.na(AnnivDates_A)))!=0) {
                # Assuming %d of LIPD as the %d of all AnnivDates produced NAs. NAs are substituted by the respective last day in month.
                nas<-which(is.na(AnnivDates_A))
                AnnivDates_B<-as.Date(timeLastDayInMonth(AnnivDates[nas]))
                AnnivDates<-sort(na.omit(append(AnnivDates_A,AnnivDates_B)))
              } else {
                AnnivDates<-AnnivDates_A
              }
            }
            AnnivDates<-sort(AnnivDates[!duplicated(AnnivDates)])
            # creating the anniversary date preceding Em
            if (AnnivDates[1]!=Em) {
              AtomVector_AD1<-as.numeric(unlist(strsplit(as.character(AnnivDates[1]),split = "-")))
              Atom1AD1<-AtomVector_AD1[1]
              Atom2AD1<-AtomVector_AD1[2]
              Atom3AD1<-AtomVector_AD1[3]
              PrevDate<-as.numeric(CppPrevDate(c(Atom1AD1,Atom2AD1,Atom3AD1,Atom1Em,Atom2Em,Atom3Em,Atom1Refer,Atom2Refer,Atom3Refer,CpY,EOM)))
              if (any(is.na(PrevDate))) {
                PrevDate<-as.Date(NA)
              } else {
                PrevDate<-as.Date(paste(PrevDate[1],PrevDate[2],PrevDate[3],sep="-"))
              }
              AnnivDates<-c(PrevDate,AnnivDates)
              AnnivDates<-sort(na.omit(AnnivDates[!duplicated(AnnivDates)]))
            }
            # creating the anniversary date succeeding Mat
            if (AnnivDates[length(AnnivDates)]!=Mat) {
              AtomVector_ADfin<-as.numeric(unlist(strsplit(as.character(AnnivDates[length(AnnivDates)]),split = "-")))
              Atom1ADfin<-AtomVector_ADfin[1]
              Atom2ADfin<-AtomVector_ADfin[2]
              Atom3ADfin<-AtomVector_ADfin[3]
              SuccDate<-as.numeric(CppSuccDate(c(Atom1ADfin,Atom2ADfin,Atom3ADfin,Atom1Mat,Atom2Mat,Atom3Mat,Atom1Refer,Atom2Refer,Atom3Refer,CpY,EOM)))
              if (any(is.na(SuccDate))) {
                SuccDate<-as.Date(NA)
              } else {
                SuccDate<-as.Date(paste(SuccDate[1],SuccDate[2],SuccDate[3],sep="-"))
              }
              AnnivDates<-c(AnnivDates,SuccDate)
              AnnivDates<-sort(na.omit(AnnivDates[!duplicated(AnnivDates)]))   # final AnnivDates created!
            }
            AD_indexes<-c(1:length(AnnivDates))-(which(AnnivDates==FIPD)-1)
            AD_List<-list(AnnivDates,AD_indexes)
            index_FIPD<-AD_List[[2]][which(AD_List[[1]]==FIPD)]
            index_LIPD<-AD_List[[2]][which(AD_List[[1]]==LIPD)]
            ###
            ### now creating the vectors RealDates and CoupDates
            ###
            #
            # AnnivDates has at least 3 elements
            # AD1 can be Em, AD_Fin can be Mat, AD_2 can be FIPD and LIPD
            #
            if (FIPD==LIPD) {   # Dealing with the case Em < FIPD = LIPD < Mat
              CoupDates<-c(FIPD,Mat)
              RealDates<-c(Em,FIPD,Mat)
            } else {            # Dealing with the case Em < FIPD < LIPD < Mat
              AD_in_RD<-as.numeric(which((!(AnnivDates<FIPD))&(!(AnnivDates>LIPD))))
              if (length(length(AD_in_RD))==0) {
                RealDates<-c(Em,Mat)
              } else {
                RealDates<-c(Em,AnnivDates[AD_in_RD],Mat)
              }
              RealDates<-RealDates[!duplicated(RealDates)]
              CoupDates<-RealDates[-1]
            }
          }
        }
        ###
        ### now using the function DIST to find E, M, f and l
        ###
        AtomVector_AD1<-as.numeric(unlist(strsplit(as.character(AnnivDates[1]),split = "-")))
        Atom1AD1<-AtomVector_AD1[1]
        Atom2AD1<-AtomVector_AD1[2]
        Atom3AD1<-AtomVector_AD1[3]
        AtomVector_AD2<-as.numeric(unlist(strsplit(as.character(AnnivDates[2]),split = "-")))
        Atom1AD2<-AtomVector_AD2[1]
        Atom2AD2<-AtomVector_AD2[2]
        Atom3AD2<-AtomVector_AD2[3]
        AtomVector_AD_PreFin<-as.numeric(unlist(strsplit(as.character(AnnivDates[length(AnnivDates)-1]),split = "-")))
        Atom1AD_PreFin<-AtomVector_AD_PreFin[1]
        Atom2AD_PreFin<-AtomVector_AD_PreFin[2]
        Atom3AD_PreFin<-AtomVector_AD_PreFin[3]
        AtomVector_AD_Fin<-as.numeric(unlist(strsplit(as.character(AnnivDates[length(AnnivDates)]),split = "-")))
        Atom1AD_Fin<-AtomVector_AD_Fin[1]
        Atom2AD_Fin<-AtomVector_AD_Fin[2]
        Atom3AD_Fin<-AtomVector_AD_Fin[3]
        if (is.element(DCC,c(1,3,5,6,8,10,11,12,15))) {
          Em_Num<-DIST(c(DCC,Atom1AD1,Atom2AD1,Atom3AD1,Atom1Em,Atom2Em,Atom3Em))[2]
          Em_Den<-DIST(c(DCC,Atom1AD1,Atom2AD1,Atom3AD1,Atom1AD2,Atom2AD2,Atom3AD2))[2]
          Mat_Num<-DIST(c(DCC,Atom1AD_PreFin,Atom2AD_PreFin,Atom3AD_PreFin,Atom1Mat,Atom2Mat,Atom3Mat))[2]
          Mat_Den<-DIST(c(DCC,Atom1AD_PreFin,Atom2AD_PreFin,Atom3AD_PreFin,Atom1AD_Fin,Atom2AD_Fin,Atom3AD_Fin))[2]
        }
        if (DCC==2|DCC==14) {
          OrigDCC<-DCC
          DCC<-2
          Em_Num<-DIST(c(DCC,Atom1AD1,Atom2AD1,Atom3AD1,Atom1AD1,Atom2AD1,Atom3AD1,Atom1AD2,Atom2AD2,Atom3AD2,
                         Atom1AD1,Atom2AD1,Atom3AD1,Atom1Em,Atom2Em,Atom3Em,Atom1AD2,Atom2AD2,Atom3AD2,
                         AD_indexes[1],AD_indexes[2],CpY))[2]
          Em_Den<-DIST(c(DCC,Atom1AD1,Atom2AD1,Atom3AD1,Atom1AD1,Atom2AD1,Atom3AD1,Atom1AD2,Atom2AD2,Atom3AD2,
                         Atom1AD2,Atom2AD2,Atom3AD2,Atom1AD2,Atom2AD2,Atom3AD2,Atom1AD2,Atom2AD2,(Atom3AD2+1),
                         AD_indexes[2],AD_indexes[2],CpY))[2]
          # Note: the 2nd line in Em_Den<-DIST(c(... should correctly be:
          #      ...,Atom1AD2,Atom2AD2,Atom3AD2,Atom1AD2,Atom2AD2,Atom3AD2,Atom1AD3,Atom2AD3,Atom3AD3,...
          # but since Em_Den calculates the distance between two coupon dates, the second summand in
          # DIST_2(PCD(t_E,AD),NCD(t_E,AD)) is 0 anyway
          Mat_Num<-DIST(c(DCC,Atom1AD_PreFin,Atom2AD_PreFin,Atom3AD_PreFin,Atom1AD_PreFin,Atom2AD_PreFin,Atom3AD_PreFin,Atom1AD_Fin,Atom2AD_Fin,Atom3AD_Fin,
                          Atom1AD_PreFin,Atom2AD_PreFin,Atom3AD_PreFin,Atom1Mat,Atom2Mat,Atom3Mat,Atom1AD_Fin,Atom2AD_Fin,Atom3AD_Fin,
                          AD_indexes[length(AD_indexes)-1],AD_indexes[length(AD_indexes)],CpY))[2]
          Mat_Den<-DIST(c(DCC,Atom1AD_PreFin,Atom2AD_PreFin,Atom3AD_PreFin,Atom1AD_PreFin,Atom2AD_PreFin,Atom3AD_PreFin,Atom1AD_Fin,Atom2AD_Fin,Atom3AD_Fin,
                          Atom1AD_Fin,Atom2AD_Fin,Atom3AD_Fin,Atom1AD_Fin,Atom2AD_Fin,Atom3AD_Fin,Atom1AD_Fin,Atom2AD_Fin,(Atom3AD_Fin+1),
                          AD_indexes[length(AD_indexes)],AD_indexes[length(AD_indexes)],CpY))[2]
          # Note: equivalently to the note above here strictly speaking one should consider
          #       the coupon anniversary date after AD_Fin in the 2nd line
          DCC<-OrigDCC
        }
        if (DCC==4) {
          Em_Num<-DIST(c(DCC,Atom1AD1,Atom2AD1,Atom3AD1,Atom1Em,Atom2Em,Atom3Em,Atom1AD2,CpY))[2]
          Em_Den<-DIST(c(DCC,Atom1AD1,Atom2AD1,Atom3AD1,Atom1AD2,Atom2AD2,Atom3AD2,Atom1AD2,CpY))[2]
          Mat_Num<-DIST(c(DCC,Atom1AD_PreFin,Atom2AD_PreFin,Atom3AD_PreFin,Atom1Mat,Atom2Mat,Atom3Mat,Atom1AD_Fin,CpY))[2]
          Mat_Den<-DIST(c(DCC,Atom1AD_PreFin,Atom2AD_PreFin,Atom3AD_PreFin,Atom1AD_Fin,Atom2AD_Fin,Atom3AD_Fin,Atom1AD_Fin,CpY))[2]
        }
        if (DCC==7) {
          Em_Num<-DIST(c(DCC,Atom1AD1,Atom2AD1,Atom3AD1,Atom1Em,Atom2Em,Atom3Em,Atom1Mat,Atom2Mat,Atom3Mat))[2]
          Em_Den<-DIST(c(DCC,Atom1AD1,Atom2AD1,Atom3AD1,Atom1AD2,Atom2AD2,Atom3AD2,Atom1Mat,Atom2Mat,Atom3Mat))[2]
          Mat_Num<-DIST(c(DCC,Atom1AD_PreFin,Atom2AD_PreFin,Atom3AD_PreFin,Atom1Mat,Atom2Mat,Atom3Mat,Atom1Mat,Atom2Mat,Atom3Mat))[2]
          Mat_Den<-DIST(c(DCC,Atom1AD_PreFin,Atom2AD_PreFin,Atom3AD_PreFin,Atom1AD_Fin,Atom2AD_Fin,Atom3AD_Fin,Atom1Mat,Atom2Mat,Atom3Mat))[2]
        }
        if (is.element(DCC,c(9,13))) {
          Em_Num<-DIST(c(DCC,Atom1AD1,Atom2AD1,Atom3AD1,Atom1Em,Atom2Em,Atom3Em,EOM))[2]
          Em_Den<-DIST(c(DCC,Atom1AD1,Atom2AD1,Atom3AD1,Atom1AD2,Atom2AD2,Atom3AD2,EOM))[2]
          Mat_Num<-DIST(c(DCC,Atom1AD_PreFin,Atom2AD_PreFin,Atom3AD_PreFin,Atom1Mat,Atom2Mat,Atom3Mat,EOM))[2]
          Mat_Den<-DIST(c(DCC,Atom1AD_PreFin,Atom2AD_PreFin,Atom3AD_PreFin,Atom1AD_Fin,Atom2AD_Fin,Atom3AD_Fin,EOM))[2]
        }
        if (DCC==16) {
          NonBus.AD1.Em<-length(which((NonBusDays.Brazil$Date>=AnnivDates[1])&(NonBusDays.Brazil$Date<Em)))
          NonBus.AD1.AD2<-length(which((NonBusDays.Brazil$Date>=AnnivDates[1])&(NonBusDays.Brazil$Date<AnnivDates[2])))
          NonBus.ADPreFin.Mat<-length(which((NonBusDays.Brazil$Date>=AnnivDates[length(AnnivDates)-1])&(NonBusDays.Brazil$Date<Mat)))
          NonBus.ADPreFin.ADFin<-length(which((NonBusDays.Brazil$Date>=AnnivDates[length(AnnivDates)-1])&(NonBusDays.Brazil$Date<AnnivDates[length(AnnivDates)])))
          Em_Num<-DIST(c(DCC,Atom1AD1,Atom2AD1,Atom3AD1,Atom1Em,Atom2Em,Atom3Em,NonBus.AD1.Em))[2]
          Em_Den<-DIST(c(DCC,Atom1AD1,Atom2AD1,Atom3AD1,Atom1AD2,Atom2AD2,Atom3AD2,NonBus.AD1.AD2))[2]
          Mat_Num<-DIST(c(DCC,Atom1AD_PreFin,Atom2AD_PreFin,Atom3AD_PreFin,Atom1Mat,Atom2Mat,Atom3Mat,NonBus.ADPreFin.Mat))[2]
          Mat_Den<-DIST(c(DCC,Atom1AD_PreFin,Atom2AD_PreFin,Atom3AD_PreFin,Atom1AD_Fin,Atom2AD_Fin,Atom3AD_Fin,NonBus.ADPreFin.ADFin))[2]
        }
        E<-Em_Num/Em_Den+AD_indexes[1]
        M<-Mat_Num/Mat_Den+AD_indexes[length(AD_indexes)-1]
        f<-1-E
        l<-M-index_LIPD
        if (index_LIPD>=1) {
          CD_indexes<-c(c(seq(1,index_LIPD,by=1)),M)
        } else {
          CD_indexes<-M
        }
        RD_indexes<-c(E,CD_indexes)
        RD_indexes<-sort(na.omit(RD_indexes[!duplicated(RD_indexes)]))
        Refer<-as.Date(paste(Atom1Refer,Atom2Refer,Atom3Refer,sep="-"))
        #
        # The following code evaluates the lengths of the first and the final coupon
        # periods f and l to determine FCPType and LCPType:
        # character string FCPType: indicates the type of the first coupon period
        # character string LCPType: indicates the type of the last coupon period
        if (round(f,3)>1) {
          FCPType<-"long"
        } else {
          if (round(f,3)<1) {
            FCPType<-"short"
          } else {
            FCPType<-"regular"
          }
        }
        if (round(l,3)>1) {
          LCPType<-"long"
        } else {
          if (round(l,3)<1) {
            LCPType<-"short"
          } else {
            LCPType<-"regular"
          }
        }
        if (Em_Mat_SameMY==1) {
          LCPType<-FCPType
          l<-f
        }
        # for DCC = {1,3,5,6,8,10,11,12,15} x is a NumericMatrix with 9 columns,
        # each row containing the vector:
        #   c(RV,Coup,DCC,Y1,M1,D1,Y2,M2,D2) with
        # Y1-M1-D1 = t_a ; Y2-M2-D2 = t_b
        #
        # for DCC = {2,14} x is a NumericMatrix with 24 columns,
        # each row containing the vector:
        #   c(RV,Coup,DCC,Y1,M1,D1,Y2,M2,D2,Y3,M3,D3,Y4,M4,D4,Y5,M5,D5,Y6,M6,D6,P,N,CpY) with
        # Y1-M1-D1 = PCD(t_a,AD) ; Y2-M2-D2 = t_a ; Y3-M3-D3 = NCD(t_a,AD)
        # Y4-M4-D4 = PCD(t_b,AD) ; Y5-M5-D5 = t_b ; Y6-M6-D6 = NCD(t_b,AD)
        # P = P(t_b,AD) ; N = N(t_a,AD)
        #
        # for DCC = 4 x is a NumericMatrix with 11 columns,
        # each row containing the vector:
        #   c(RV,Coup,DCC,Y1,M1,D1,Y2,M2,D2,Y3,CpY) with
        # Y1-M1-D1 = t_a ; Y2-M2-D2 = t_b ; Y3-M3-D3 = t_c
        #
        # for DCC = 7 x is a NumericMatrix with 12 columns,
        # each row containing the vector:
        #   c(RV,Coup,DCC,Y1,M1,D1,Y2,M2,D2,Y3,M3,D3) with
        # Y1-M1-D1 = t_a ; Y2-M2-D2 = t_b ; Y3-M3-D3 = t_M
        #
        # for DCC = {9,13} x is a NumericMatrix with 10 columns,
        # each row containing the vector:
        #   c(RV,Coup,DCC,Y1,M1,D1,Y2,M2,D2,EOM) with
        # Y1-M1-D1 = t_a ; Y2-M2-D2 = t_b
        #
        # for DCC = 16 x is a NumericMatrix with 10 columns,
        # each row containing the vector:
        #   c(RV,Coup,DCC,Y1,M1,D1,Y2,M2,D2,NonBus) with
        # Y1-M1-D1 = t_a ; Y2-M2-D2 = t_b
        # NonBus = # of non-business days between t_a (incl) and t_b (excl)
        #
        if (!(is.na(Coup))) {
          if (is.element(DCC,c(1,3,5,6,8,10,11,12,15))) {
            RVVec<-rep(RV,length(RealDates)-1)
            CoupVec<-rep((Coup/100),length(RealDates)-1)
            DCCVec<-rep(DCC,length(RealDates)-1)
            D1<-RealDates[-length(RealDates)]
            D2<-RealDates[-1]
            D1Atoms<-as.numeric(unlist(strsplit(as.character(D1),split="-")))
            D2Atoms<-as.numeric(unlist(strsplit(as.character(D2),split="-")))
            D1Matrix<-matrix(t(D1Atoms),nrow=length(RealDates)-1,ncol=3,byrow=TRUE)
            D2Matrix<-matrix(t(D2Atoms),nrow=length(RealDates)-1,ncol=3,byrow=TRUE)
            CoupSchedMatrix<-cbind(RVVec,CoupVec,DCCVec,D1Matrix,D2Matrix)
            CoupPayments<-PayCalc(CoupSchedMatrix)
          }
          if (DCC==16) {
            NonBus.Em.CD1<-length(which((NonBusDays.Brazil$Date>=Em)&(NonBusDays.Brazil$Date<CoupDates[1])))
            AtomVector_CD1<-as.numeric(unlist(strsplit(as.character(CoupDates[1]),split = "-")))
            Atom1CD1<-AtomVector_CD1[1]
            Atom2CD1<-AtomVector_CD1[2]
            Atom3CD1<-AtomVector_CD1[3]
            Exp.FirstCoup<-DIST(c(DCC,Atom1Em,Atom2Em,Atom3Em,Atom1CD1,Atom2CD1,Atom3CD1,NonBus.Em.CD1))[2]
            FirstCoup<-RV*(((1+Coup/100)^(Exp.FirstCoup))-1)
            if (length(RealDates)==2) {
              CoupPayments<-FirstCoup
            } else {
              RegCoup<-RV*(((1+Coup/100)^(1/CpY))-1)
              NonBus.CD_PreFin.Mat<-length(which((NonBusDays.Brazil$Date>=CoupDates[length(CoupDates)-1])&(NonBusDays.Brazil$Date<Mat)))
              AtomVector_CD_PreFin<-as.numeric(unlist(strsplit(as.character(CoupDates[length(CoupDates)-1]),split = "-")))
              Atom1CD_PreFin<-AtomVector_CD_PreFin[1]
              Atom2CD_PreFin<-AtomVector_CD_PreFin[2]
              Atom3CD_PreFin<-AtomVector_CD_PreFin[3]
              Exp.FinalCoup<-DIST(c(DCC,Atom1CD_PreFin,Atom2CD_PreFin,Atom3CD_PreFin,Atom1Mat,Atom2Mat,Atom3Mat,NonBus.CD_PreFin.Mat))[2]
              FinalCoup<-RV*(((1+Coup/100)^(Exp.FinalCoup))-1)
              CoupPayments<-c(FirstCoup,rep(RegCoup,(length(CoupDates)-2)),FinalCoup)
            }
          }
          if (DCC==2|DCC==14) {
            # creating the anniversary date preceding AD1
            AtomVector_AD1<-as.numeric(unlist(strsplit(as.character(AnnivDates[1]),split = "-")))
            Atom1AD1<-AtomVector_AD1[1]
            Atom2AD1<-AtomVector_AD1[2]
            Atom3AD1<-AtomVector_AD1[3]
            PrevDate<-as.numeric(CppPrevDate(c(Atom1AD1,Atom2AD1,Atom3AD1,Atom1AD1,Atom2AD1,Atom3AD1,Atom1Refer,Atom2Refer,Atom3Refer,CpY,EOM)))
            PrevDate<-as.Date(paste(PrevDate[1],PrevDate[2],PrevDate[3],sep="-"))
            AnnivDates_PaySched<-c(PrevDate,AnnivDates)
            AnnivDates_PaySched<-sort(na.omit(AnnivDates_PaySched[!duplicated(AnnivDates_PaySched)]))
            # creating the anniversary date succeeding ADfin
            AtomVector_ADfin<-as.numeric(unlist(strsplit(as.character(AnnivDates[length(AnnivDates)]),split = "-")))
            Atom1ADfin<-AtomVector_ADfin[1]
            Atom2ADfin<-AtomVector_ADfin[2]
            Atom3ADfin<-AtomVector_ADfin[3]
            SuccDate<-as.numeric(CppSuccDate(c(Atom1ADfin,Atom2ADfin,Atom3ADfin,Atom1ADfin,Atom2ADfin,Atom3ADfin,Atom1Refer,Atom2Refer,Atom3Refer,CpY,EOM)))
            SuccDate<-as.Date(paste(SuccDate[1],SuccDate[2],SuccDate[3],sep="-"))
            AnnivDates_PaySched<-c(AnnivDates_PaySched,SuccDate)
            AnnivDates_PaySched<-sort(na.omit(AnnivDates_PaySched[!duplicated(AnnivDates_PaySched)]))
            #   c(RV,Coup,DCC,Y1,M1,D1,Y2,M2,D2,Y3,M3,D3,Y4,M4,D4,Y5,M5,D5,Y6,M6,D6,P,N,CpY)
            RVVec<-rep(RV,length(RealDates)-1)
            CoupVec<-rep((Coup/100),length(RealDates)-1)
            DCCVec<-rep(DCC,length(RealDates)-1)
            D1<-RealDates[-length(RealDates)]
            PrevDates_D1<-as.Date(unlist(lapply(D1,PCD,AnnivDates_PaySched)),DateOrigin)
            SuccDates_D1<-as.Date(unlist(lapply(D1,NCD,AnnivDates_PaySched)),DateOrigin)
            D2<-RealDates[-1]
            PrevDates_D2<-as.Date(unlist(lapply(D2,PCD,AnnivDates_PaySched)),DateOrigin)
            SuccDates_D2<-as.Date(unlist(lapply(D2,NCD,AnnivDates_PaySched)),DateOrigin)
            P_D2<-rep(AD_List[[2]],length(PrevDates_D2))[which(unlist(lapply(PrevDates_D2, `==`, AD_List[[1]]))==TRUE)]
            N_D1<-rep(AD_List[[2]],length(SuccDates_D1))[which(unlist(lapply(SuccDates_D1, `==`, AD_List[[1]]))==TRUE)]
            CpYVec<-rep(CpY,length(RealDates)-1)
            # extracting atoms from dates
            PrevDates_D1_Atoms<-as.numeric(unlist(strsplit(as.character(PrevDates_D1),split="-")))
            D1Atoms<-as.numeric(unlist(strsplit(as.character(D1),split="-")))
            SuccDates_D1_Atoms<-as.numeric(unlist(strsplit(as.character(SuccDates_D1),split="-")))
            PrevDates_D2_Atoms<-as.numeric(unlist(strsplit(as.character(PrevDates_D2),split="-")))
            D2Atoms<-as.numeric(unlist(strsplit(as.character(D2),split="-")))
            SuccDates_D2_Atoms<-as.numeric(unlist(strsplit(as.character(SuccDates_D2),split="-")))
            # constructing matrices of the date atoms
            PrevDates_D1_Matrix<-matrix(t(PrevDates_D1_Atoms),nrow=length(RealDates)-1,ncol=3,byrow=TRUE)
            D1Matrix<-matrix(t(D1Atoms),nrow=length(RealDates)-1,ncol=3,byrow=TRUE)
            SuccDates_D1_Matrix<-matrix(t(SuccDates_D1_Atoms),nrow=length(RealDates)-1,ncol=3,byrow=TRUE)
            PrevDates_D2_Matrix<-matrix(t(PrevDates_D2_Atoms),nrow=length(RealDates)-1,ncol=3,byrow=TRUE)
            D2Matrix<-matrix(t(D2Atoms),nrow=length(RealDates)-1,ncol=3,byrow=TRUE)
            SuccDates_D2_Matrix<-matrix(t(SuccDates_D2_Atoms),nrow=length(RealDates)-1,ncol=3,byrow=TRUE)
            # binding everything into one matrix
            CoupSchedMatrix<-cbind(RVVec,CoupVec,DCCVec,PrevDates_D1_Matrix,D1Matrix,SuccDates_D1_Matrix,
                                   PrevDates_D2_Matrix,D2Matrix,SuccDates_D2_Matrix,P_D2,N_D1,CpYVec)
            CoupPayments<-PayCalc(CoupSchedMatrix)
          }
          if (DCC==4) {
            # c(RV,Coup,DCC,Y1,M1,D1,Y2,M2,D2,Y3,CpY)
            RVVec<-rep(RV,length(RealDates)-1)
            CoupVec<-rep((Coup/100),length(RealDates)-1)
            DCCVec<-rep(DCC,length(RealDates)-1)
            D1<-RealDates[-length(RealDates)]
            D2<-RealDates[-1]
            D1Atoms<-as.numeric(unlist(strsplit(as.character(D1),split="-")))
            D2Atoms<-as.numeric(unlist(strsplit(as.character(D2),split="-")))
            D1Matrix<-matrix(t(D1Atoms),nrow=length(RealDates)-1,ncol=3,byrow=TRUE)
            D2Matrix<-matrix(t(D2Atoms),nrow=length(RealDates)-1,ncol=3,byrow=TRUE)
            Y3Vector<-D2Matrix[,1]
            CpYVec<-rep(CpY,length(RealDates)-1)
            CoupSchedMatrix<-cbind(RVVec,CoupVec,DCCVec,D1Matrix,D2Matrix,Y3Vector,CpYVec)
            CoupPayments<-PayCalc(CoupSchedMatrix)
          }
          if (DCC==7) {
            # c(RV,Coup,DCC,Y1,M1,D1,Y2,M2,D2,Y3,M3,D3)
            RVVec<-rep(RV,length(RealDates)-1)
            CoupVec<-rep((Coup/100),length(RealDates)-1)
            DCCVec<-rep(DCC,length(RealDates)-1)
            D1<-RealDates[-length(RealDates)]
            D2<-RealDates[-1]
            MatVec<-rep(Mat,length(RealDates)-1)
            D1Atoms<-as.numeric(unlist(strsplit(as.character(D1),split="-")))
            D2Atoms<-as.numeric(unlist(strsplit(as.character(D2),split="-")))
            MatAtoms<-as.numeric(unlist(strsplit(as.character(MatVec),split="-")))
            D1Matrix<-matrix(t(D1Atoms),nrow=length(RealDates)-1,ncol=3,byrow=TRUE)
            D2Matrix<-matrix(t(D2Atoms),nrow=length(RealDates)-1,ncol=3,byrow=TRUE)
            MatMatrix<-matrix(t(MatAtoms),nrow=length(RealDates)-1,ncol=3,byrow=TRUE)
            CoupSchedMatrix<-cbind(RVVec,CoupVec,DCCVec,D1Matrix,D2Matrix,MatMatrix)
            CoupPayments<-PayCalc(CoupSchedMatrix)
          }
          if (is.element(DCC,c(9,13))) {
            #   c(RV,Coup,DCC,Y1,M1,D1,Y2,M2,D2,EOM) with
            RVVec<-rep(RV,length(RealDates)-1)
            CoupVec<-rep((Coup/100),length(RealDates)-1)
            DCCVec<-rep(DCC,length(RealDates)-1)
            D1<-RealDates[-length(RealDates)]
            D2<-RealDates[-1]
            EOMVec<-rep(EOM,length(RealDates)-1)
            D1Atoms<-as.numeric(unlist(strsplit(as.character(D1),split="-")))
            D2Atoms<-as.numeric(unlist(strsplit(as.character(D2),split="-")))
            EOMAtoms<-as.numeric(unlist(strsplit(as.character(EOMVec),split="-")))
            D1Matrix<-matrix(t(D1Atoms),nrow=length(RealDates)-1,ncol=3,byrow=TRUE)
            D2Matrix<-matrix(t(D2Atoms),nrow=length(RealDates)-1,ncol=3,byrow=TRUE)
            EOMMatrix<-matrix(t(EOMAtoms),nrow=length(RealDates)-1,ncol=3,byrow=TRUE)
            CoupSchedMatrix<-cbind(RVVec,CoupVec,DCCVec,D1Matrix,D2Matrix,EOMMatrix)
            CoupPayments<-PayCalc(CoupSchedMatrix)
          }
          if ((!(RegCF.equal==0))&(!(DCC==16))) {
            if (length(CoupPayments)>2) {
              CoupPayments<-
                c(CoupPayments[1],rep(RV*(Coup/(100*CpY)),(length(CoupPayments)-2)),CoupPayments[length(CoupPayments)])
              PaySched<-data.frame(CoupDates,CoupPayments)
            } else {
              PaySched<-data.frame(CoupDates,CoupPayments)
            }
          } else {
            PaySched<-data.frame(CoupDates,CoupPayments)
          }
        }
      }
    }
  }
  Warnings<-data.frame(Em_FIAD_differ=Em_FIAD_differ,EmMatMissing=EmMatMissing,CpYOverride=CpYOverride,RV_set100percent=RV_set100percent,
                       NegLifeFlag=NegLifeFlag,ZeroFlag=ZeroFlag,Em_Mat_SameMY=Em_Mat_SameMY,ChronErrorFlag=ChronErrorFlag,
                       FIPD_LIPD_equal=FIPD_LIPD_equal,IPD_CpY_Corrupt=IPD_CpY_Corrupt,EOM_Deviation=EOM_Deviation,EOMOverride=EOMOverride,
                       DCCOverride=DCCOverride,NoCoups=NoCoups)
  Traits<-data.frame(DateOrigin=DateOrigin,CpY=CpY,FIAD=FIAD,Em=Em,Em_Orig=Em_Orig,FIPD=FIPD,FIPD_Orig=FIPD_Orig,est_FIPD=est_FIPD,LIPD=LIPD,
                     LIPD_Orig=LIPD_Orig,est_LIPD=est_LIPD,Mat=Mat,Refer=Refer,FCPType=FCPType,FCPLength=f,LCPType=LCPType,LCPLength=l,Par=RV,
                     CouponInPercent.p.a=Coup,DayCountConvention=DCC,EOM_Orig=EOM_Orig,est_EOM=est_EOM,EOM_used=EOM)
  RealDates<-append(RealDates,rep(NA,max(length(RealDates),length(CoupDates),length(AnnivDates))-length(RealDates)))
  RD_indexes<-append(RD_indexes,rep(NA,max(length(RealDates),length(CoupDates),length(AnnivDates))-length(RD_indexes)))
  CoupDates<-append(CoupDates,rep(NA,max(length(RealDates),length(CoupDates),length(AnnivDates))-length(CoupDates)))
  CD_indexes<-append(CD_indexes,rep(NA,max(length(RealDates),length(CoupDates),length(AnnivDates))-length(CD_indexes)))
  AnnivDates<-append(AnnivDates,rep(NA,max(length(RealDates),length(CoupDates),length(AnnivDates))-length(AnnivDates)))
  AD_indexes<-append(AD_indexes,rep(NA,max(length(RealDates),length(CoupDates),length(AnnivDates))-length(AnnivDates)))
  DateVectors<-data.frame(RealDates=RealDates,RD_indexes=RD_indexes,CoupDates=CoupDates,CD_indexes=CD_indexes,AnnivDates=AnnivDates,AD_indexes=AD_indexes)
  if (all(is.na(PaySched))) {
    AnnivDates_Out<-list(Warnings=Warnings,Traits=Traits,DateVectors=DateVectors)
  } else {
    AnnivDates_Out<-list(Warnings=Warnings,Traits=Traits,DateVectors=DateVectors,PaySched=PaySched)
  }
  return(AnnivDates_Out)
}

