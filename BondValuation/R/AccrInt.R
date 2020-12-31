#' AccrInt (calculation of accrued interest)
#'
#' \bold{AccrInt} returns the amount of interest accrued from some starting date
#' up to some end date and the number of days of interest on the end date.
#'
#'     \tabular{cl}{
#'       \code{DCC} \tab required input \cr
#'           _____________________ \tab _____________________________________________ \cr
#'       1,3,5,6,8,10,11,12,15,16 \tab \code{StartDate}, \code{EndDate},
#'                                     \code{Coup}, \code{DCC}, \code{RV} \cr
#'                           2,14 \tab \code{StartDate}, \code{EndDate},
#'                                     \code{Coup}, \code{DCC}, \code{RV},
#'                                     \code{CpY}, \code{EOM} \cr
#'                              4 \tab \code{StartDate}, \code{EndDate},
#'                                     \code{Coup}, \code{DCC}, \code{RV},
#'                                     \code{CpY}, \code{EOM},
#'                                     \code{YearNCP} \cr
#'                              7 \tab \code{StartDate}, \code{EndDate},
#'                                     \code{Coup}, \code{DCC}, \code{RV},
#'                                     \code{Mat} \cr
#'                           9,13 \tab \code{StartDate}, \code{EndDate},
#'                                     \code{Coup}, \code{DCC}, \code{RV},
#'                                     \code{EOM} \cr
#'            =================== \tab ======================================== \cr
#'     }
#'
#' Assuming that there is no accrued interest on \code{StartDate} the function
#' \bold{AccrInt} computes the amount of interest accrued up to \code{EndDate}
#' under the terms of the specified day count convention \code{DCC}. The function
#' returns a list of two numerics \code{AccrInt}, and \code{DaysAccrued}.
#' If \code{InputCheck = 1} the input variables are checked for the correct
#' format. The core feature of this function is the proper handling of the
#' \emph{\bold{day count conventions}} presented below. The type of the day
#' count convention determines the amount of the accrued interest that has
#' to be paid by the buyer in the secondary market if the settlement
#' takes place between two coupon payment dates.
#'
#' \itemize{
#'   \item Many different day count conventions are used in the market.
#'         Since there is no central authority that develops these conventions
#'         there is no standardized nomenclature. The tables below provide
#'         alternative names that often are used for the respective conventions.
#'         Type \code{View(List.DCC)} for a list of the day count methods
#'         currently implemented.
#'   \item Detailed descriptions of the conventions and their application may
#'         be found in Djatschenko (2018), and the other provided references.
#' }
#'
#' \bold{Day Count Conventions}
#'
#' \describe{
#'   \item{}{
#'     \tabular{cccl}{
#'                    \tab   \tab     \tab \bold{Actual/Actual (ISDA)}                      \cr
#'         ___________\tab | \tab ___ \tab ________________________________________________ \cr
#'                DCC \tab | \tab  =  \tab 1                                                \cr
#'         ___________\tab | \tab ___ \tab ________________________________________________ \cr
#'        other names \tab | \tab     \tab Actual/Actual, Act/Act, Act/Act (ISDA)           \cr
#'         ___________\tab | \tab ___ \tab ________________________________________________ \cr
#'         references \tab | \tab     \tab ISDA (1998); ISDA (2006) section 4.16 (b)        \cr
#'          ==========\tab | \tab === \tab ===========================================      \cr
#'     }
#'   }
#'
#'   \item{}{
#'     \tabular{cccl}{
#'                    \tab   \tab     \tab \bold{Actual/Actual (ICMA)}                      \cr
#'         ___________\tab | \tab ___ \tab ________________________________________________ \cr
#'                DCC \tab | \tab  =  \tab 2                                                \cr
#'         ___________\tab | \tab ___ \tab ________________________________________________ \cr
#'        other names \tab | \tab     \tab Actual/Actual (ISMA), Act/Act (ISMA),            \cr
#'                    \tab | \tab     \tab Act/Act (ICMA), ISMA-99                          \cr
#'         ___________\tab | \tab ___ \tab ________________________________________________ \cr
#'         references \tab | \tab     \tab ICMA Rule 251; ISDA (2006) section 4.16 (c);     \cr
#'                    \tab | \tab     \tab SWX (2003)                                       \cr
#'          ==========\tab | \tab === \tab ===========================================      \cr
#'     }
#'   }
#'
#'   \item{}{
#'     \tabular{cccl}{
#'                    \tab   \tab     \tab \bold{Actual/Actual (AFB)}                       \cr
#'         ___________\tab | \tab ___ \tab ________________________________________________ \cr
#'                DCC \tab | \tab  =  \tab 3                                                \cr
#'         ___________\tab | \tab ___ \tab ________________________________________________ \cr
#'        other names \tab | \tab     \tab AFB Method, Actual/Actual (Euro),                \cr
#'                    \tab | \tab     \tab Actual/Actual AFB FBF, ACT/365-366 (leap day)    \cr
#'         ___________\tab | \tab ___ \tab ________________________________________________ \cr
#'         references \tab | \tab     \tab ISDA (1998); EBF (2004)                          \cr
#'          ==========\tab | \tab === \tab ===========================================      \cr
#'     }
#'   }
#'
#'   \item{}{
#'     \tabular{cccl}{
#'                    \tab   \tab     \tab \bold{Actual/365L}                               \cr
#'         ___________\tab | \tab ___ \tab ________________________________________________ \cr
#'                DCC \tab | \tab  =  \tab 4                                                \cr
#'         ___________\tab | \tab ___ \tab ________________________________________________ \cr
#'        other names \tab | \tab     \tab Act/365-366, ISMA-Year                           \cr
#'         ___________\tab | \tab ___ \tab ________________________________________________ \cr
#'         references \tab | \tab     \tab ICMA Rule 251; SWX (2003)                        \cr
#'          ==========\tab | \tab === \tab ===========================================      \cr
#'     }
#'   }
#'
#'   \item{}{
#'     \tabular{cccl}{
#'                    \tab   \tab     \tab \bold{30/360}                                    \cr
#'         ___________\tab | \tab ___ \tab ________________________________________________ \cr
#'                DCC \tab | \tab  =  \tab 5                                                \cr
#'         ___________\tab | \tab ___ \tab ________________________________________________ \cr
#'        other names \tab | \tab     \tab 360/360, Bond Basis, 30/360 ISDA                 \cr
#'         ___________\tab | \tab ___ \tab ________________________________________________ \cr
#'         references \tab | \tab     \tab ISDA (2006) section 4.16 (f);                    \cr
#'                    \tab | \tab     \tab MSRB (2017) Rule G-33                            \cr
#'          ==========\tab | \tab === \tab ===========================================      \cr
#'     }
#'   }
#'
#'   \item{}{
#'     \tabular{cccl}{
#'                    \tab   \tab     \tab \bold{30E/360}                                   \cr
#'         ___________\tab | \tab ___ \tab ________________________________________________ \cr
#'                DCC \tab | \tab  =  \tab 6                                                \cr
#'         ___________\tab | \tab ___ \tab ________________________________________________ \cr
#'        other names \tab | \tab     \tab Eurobond Basis, Special German (30S/360),        \cr
#'                    \tab | \tab     \tab ISMA-30/360                                      \cr
#'         ___________\tab | \tab ___ \tab ________________________________________________ \cr
#'         references \tab | \tab     \tab ICMA Rule 251; ISDA (2006) section 4.16 (g);     \cr
#'                    \tab | \tab     \tab SWX (2003)                                       \cr
#'          ==========\tab | \tab === \tab ===========================================      \cr
#'     }
#'   }
#'
#'   \item{}{
#'     \tabular{cccl}{
#'                    \tab   \tab     \tab \bold{30E/360 (ISDA)}                            \cr
#'         ___________\tab | \tab ___ \tab ________________________________________________ \cr
#'                DCC \tab | \tab  =  \tab 7                                                \cr
#'         ___________\tab | \tab ___ \tab ________________________________________________ \cr
#'        other names \tab | \tab     \tab none                                             \cr
#'         ___________\tab | \tab ___ \tab ________________________________________________ \cr
#'         references \tab | \tab     \tab ISDA (2006) section 4.16 (h)                     \cr
#'          ==========\tab | \tab === \tab ===========================================      \cr
#'     }
#'   }
#'
#'   \item{}{
#'     \tabular{cccl}{
#'                    \tab   \tab     \tab \bold{30/360 (German)}                           \cr
#'         ___________\tab | \tab ___ \tab ________________________________________________ \cr
#'                DCC \tab | \tab  =  \tab 8                                                \cr
#'         ___________\tab | \tab ___ \tab ________________________________________________ \cr
#'        other names \tab | \tab     \tab 360/360 (German Master); German (30/360)         \cr
#'         ___________\tab | \tab ___ \tab ________________________________________________ \cr
#'         references \tab | \tab     \tab EBF (2004); SWX (2003)                           \cr
#'          ==========\tab | \tab === \tab ===========================================      \cr
#'     }
#'   }
#'
#'   \item{}{
#'     \tabular{cccl}{
#'                    \tab   \tab     \tab \bold{30/360 US}                                 \cr
#'         ___________\tab | \tab ___ \tab ________________________________________________ \cr
#'                DCC \tab | \tab  =  \tab 9                                                \cr
#'         ___________\tab | \tab ___ \tab ________________________________________________ \cr
#'        other names \tab | \tab     \tab 30/360, US (30U/360), 30/360 (SIA)               \cr
#'         ___________\tab | \tab ___ \tab ________________________________________________ \cr
#'         references \tab | \tab     \tab Mayle (1993); SWX (2003)                         \cr
#'          ==========\tab | \tab === \tab ===========================================      \cr
#'     }
#'   }
#'
#'   \item{}{
#'     \tabular{cccl}{
#'                    \tab   \tab     \tab \bold{Actual/365 (Fixed)}                        \cr
#'         ___________\tab | \tab ___ \tab ________________________________________________ \cr
#'                DCC \tab | \tab  =  \tab 10                                               \cr
#'         ___________\tab | \tab ___ \tab ________________________________________________ \cr
#'        other names \tab | \tab     \tab Act/365 (Fixed), A/365 (Fixed), A/365F, English  \cr
#'         ___________\tab | \tab ___ \tab ________________________________________________ \cr
#'         references \tab | \tab     \tab ISDA (2006) section 4.16 (d); SWX (2003)         \cr
#'          ==========\tab | \tab === \tab ===========================================      \cr
#'     }
#'   }
#'
#'   \item{}{
#'     \tabular{cccl}{
#'                    \tab   \tab     \tab \bold{Actual(NL)/365}                            \cr
#'         ___________\tab | \tab ___ \tab ________________________________________________ \cr
#'                DCC \tab | \tab  =  \tab 11                                               \cr
#'         ___________\tab | \tab ___ \tab ________________________________________________ \cr
#'        other names \tab | \tab     \tab Act(No Leap Year)/365                            \cr
#'         ___________\tab | \tab ___ \tab ________________________________________________ \cr
#'         references \tab | \tab     \tab Krgin (2002); Thomson Reuters EIKON              \cr
#'          ==========\tab | \tab === \tab ===========================================      \cr
#'                    \tab   \tab     \tab                                                  \cr
#'                    \tab   \tab     \tab                                                  \cr
#'     }
#'   }
#'
#'   \item{}{
#'     \tabular{cccl}{
#'                    \tab   \tab     \tab \bold{Actual/360}                                \cr
#'         ___________\tab | \tab ___ \tab ________________________________________________ \cr
#'                DCC \tab | \tab  =  \tab 12                                               \cr
#'         ___________\tab | \tab ___ \tab ________________________________________________ \cr
#'        other names \tab | \tab     \tab Act/360, A/360, French                           \cr
#'         ___________\tab | \tab ___ \tab ________________________________________________ \cr
#'         references \tab | \tab     \tab ISDA (2006) section 4.16 (e); SWX (2003)         \cr
#'          ==========\tab | \tab === \tab ===========================================      \cr
#'     }
#'   }
#'
#'   \item{}{
#'     \tabular{cccl}{
#'                    \tab   \tab     \tab \bold{30/365}                                    \cr
#'         ___________\tab | \tab ___ \tab ________________________________________________ \cr
#'                DCC \tab | \tab  =  \tab 13                                               \cr
#'         ___________\tab | \tab ___ \tab ________________________________________________ \cr
#'         references \tab | \tab     \tab Krgin (2002); Thomson Reuters EIKON              \cr
#'          ==========\tab | \tab === \tab ===========================================      \cr
#'     }
#'   }
#'
#'   \item{}{
#'     \tabular{cccl}{
#'                    \tab   \tab     \tab \bold{Act/365 (Canadian Bond)}                   \cr
#'         ___________\tab | \tab ___ \tab ________________________________________________ \cr
#'                DCC \tab | \tab  =  \tab 14                                               \cr
#'         ___________\tab | \tab ___ \tab ________________________________________________ \cr
#'         references \tab | \tab     \tab IIAC (2018); Thomson Reuters EIKON               \cr
#'          ==========\tab | \tab === \tab ===========================================      \cr
#'     }
#'   }
#'
#'   \item{}{
#'     \tabular{cccl}{
#'                    \tab   \tab     \tab \bold{Act/364}                                   \cr
#'         ___________\tab | \tab ___ \tab ________________________________________________ \cr
#'                DCC \tab | \tab  =  \tab 15                                               \cr
#'         ___________\tab | \tab ___ \tab ________________________________________________ \cr
#'         references \tab | \tab     \tab Thomson Reuters EIKON                            \cr
#'          ==========\tab | \tab === \tab ===========================================      \cr
#'     }
#'   }
#'
#'   \item{}{
#'     \tabular{cccl}{
#'                    \tab   \tab     \tab \bold{BusDay/252 (Brazilian)}                    \cr
#'         ___________\tab | \tab ___ \tab ________________________________________________ \cr
#'                DCC \tab | \tab  =  \tab 16                                               \cr
#'         ___________\tab | \tab ___ \tab ________________________________________________ \cr
#'        other names \tab | \tab     \tab BUS/252, BD/252                                  \cr
#'         ___________\tab | \tab ___ \tab ________________________________________________ \cr
#'         references \tab | \tab     \tab Caputo Silva et al. (2010),                      \cr
#'                    \tab | \tab     \tab Itau Unibanco S.A. (2017)                        \cr
#'          ==========\tab | \tab === \tab ===========================================      \cr
#'     }
#'   }
#'
#' }
#'
#' @param StartDate Calendar date on which interest accrual starts. Date class object with format "\%Y-\%m-\%d". (required)
#' @param EndDate Calendar date up to which interest accrues. Date class object with format "\%Y-\%m-\%d". (required)
#' @param Coup Nominal interest rate per year in percent. (required)
#' @param DCC The day count convention for interest accrual. (required)
#' @param RV The redemption value of the bond. Default: 100.
#' @param CpY Number of interest payments per year (non-negative integer; element of the set
#'        \{1,2,3,4,6,12\}. Default: 2.
#' @param Mat So-called "maturity date" i.e. date on which the redemption value and the final interest
#'        are paid. Date class object with format "\%Y-\%m-\%d".
#' @param YearNCP Year figure of the next coupon payment date after \code{EndDate}.
#' @param EOM Boolean indicating whether the bond follows the End-of-Month rule.
#' @param DateOrigin Determines the starting point for the daycount in "Date" objects.
#'        Default: "1970-01-01".
#' @param InputCheck If 1, the input variables are checked for the correct format. Default: 1.
#'
#' @return
#'   \describe{
#'     \item{AccrInt}{
#'       Accrued interest on \code{EndDate}, given the other characteristics.
#'     }
#'     \item{DaysAccrued}{
#'       The number of days of interest from \code{StartDate} to \code{EndDate}.
#'     }
#'   }
#'
#' @references
#' \enumerate{
#'   \item{Banking Federation of the European Union (EBF), 2004, Master Agreement for Financial Transactions - Supplement to the Derivatives Annex - Interest Rate Transactions.}
#'   \item{Caputo Silva, Anderson, Lena Oliveira de Carvalho, and Octavio Ladeira de Medeiros, 2010, \emph{Public Debt: The Brazilian Experience} (National Treasury Secretariat and World Bank, Brasilia, BR).}
#'   \item{Djatschenko, Wadim, The Nitty Gritty of Bond Valuation: A Generalized Methodology for Fixed Coupon Bond Analysis Allowing for Irregular Periods and Various Day Count Conventions (November 5, 2018). Available at SSRN: https://ssrn.com/abstract=3205167.}
#'   \item{International Capital Market Association (ICMA), 2010, Rule 251 Accrued Interest Calculation - Excerpt from ICMA's Rules and Recommendations.}
#'   \item{Investment Industry Association of Canada (IIAC), 2018, Canadian Conventions in Fixed Income Markets - A Reference Document of Fixed Income Securities Formulas and Practices; Release: 1.3.}
#'   \item{International Swaps and Derivatives Association (ISDA), Inc., 1998, "EMU and Market Conventions: Recent Developments".}
#'   \item{International Swaps and Derivatives Association (ISDA), 2006, Inc., \emph{2006 ISDA Definitions.}, New York.}
#'   \item{Itau Unibanco S.A., 2017, Brazilian Sovereign Fixed Income and Foreign Exchange Markets - Handbook (First Edition).}
#'   \item{Krgin, Dragomir, 2002, The Handbook of Global Fixed Income Calculations. (Wiley, New York).}
#'   \item{Mayle, Jan, 1993, Standard Securities Calculation Methods: Fixed Income Securities Formulas for Price, Yield, and Accrued Interest, volume 1, New York: Securities Industry Association, third edition.}
#'   \item{Municipal Securities Rulemaking Board (MSRB), 2017, MSRB Rule Book, Washington, DC: Municipal Securities Rulemaking Board.}
#'   \item{SWX Swiss Exchange and D. Christie, 2003, "Accrued Interest & Yield Calculations and Determination of Holiday Calendars".}
#' }
#'
#' @examples
#' StartDate<-rep(as.Date("2011-08-31"),16)
#' EndDate<-rep(as.Date("2012-02-29"),16)
#' Coup<-rep(5.25,16)
#' DCC<-seq(1,16)
#' RV<-rep(10000,16)
#' CpY<-rep(2,16)
#' Mat<-rep(as.Date("2021-08-31"),16)
#' YearNCP<-rep(2012,16)
#' EOM<-rep(1,16)
#'
#' DCC_Comparison<-data.frame(StartDate,EndDate,Coup,DCC,RV,CpY,Mat,YearNCP,EOM)
#'
#' AccrIntOutput<-apply(DCC_Comparison[,c('StartDate','EndDate','Coup','DCC',
#' 'RV','CpY','Mat','YearNCP','EOM')],1,function(y) AccrInt(y[1],y[2],y[3],
#' y[4],y[5],y[6],y[7],y[8],y[9]))
#' # warnings are due to apply's conversion of the variables' classes in
#' # DCC_Comparison to class "character"
#' Accrued_Interest<-do.call(rbind,lapply(AccrIntOutput, function(x) x[[1]]))
#' Days_Accrued<-do.call(rbind,lapply(AccrIntOutput, function(x) x[[2]]))
#' DCC_Comparison<-cbind(DCC_Comparison,Accrued_Interest,Days_Accrued)
#' DCC_Comparison
#'
#'
#' @export
AccrInt<-function(StartDate=as.Date(NA),EndDate=as.Date(NA),Coup=as.numeric(NA),DCC=as.numeric(NA),RV=as.numeric(NA),CpY=as.numeric(NA),Mat=as.Date(NA),YearNCP=as.Date(NA),EOM=as.numeric(NA),DateOrigin=as.Date("1970-01-01"),InputCheck=1) {
  NAccr<-as.numeric(NA)
  AccrInt<-as.numeric(NA)
  if (InputCheck==1) {
    CheckedInput<-InputFormatCheck(StartDate=StartDate,EndDate=EndDate,Coup=Coup,DCC=DCC,RV=RV,CpY=CpY,Mat=Mat,YearNCP=YearNCP,EOM=EOM,DateOrigin=DateOrigin)
    StartDate<-CheckedInput$StartDate
    EndDate<-CheckedInput$EndDate
    Coup<-CheckedInput$Coup
    DCC<-CheckedInput$DCC
    RV<-CheckedInput$RV
    CpY<-CheckedInput$CpY
    Mat<-CheckedInput$Mat
    YearNCP<-CheckedInput$YearNCP
    EOM<-CheckedInput$EOM
    DateOrigin<-CheckedInput$DateOrigin
  }
  if ((missing(StartDate))|(is.na(StartDate))) {
    warning("The supplied StartDate is NA or cannot be processed. NA created!")
  } else {
    if ((missing(EndDate))|(is.na(EndDate))) {
      warning("The supplied EndDate is NA or cannot be processed. NA created!")
    } else {
      if (EndDate<StartDate) {
        warning("The supplied EndDate is prior to the supplied StartDate. NA created!")
      } else {
        if ((missing(Coup))|(is.na(Coup))) {
          warning("The supplied interest rate p.a. (Coup) is NA or cannot be processed. NA created!")
        } else {
          # If DCC is not provided or NA or not element of {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16}, the following code sets it 2 (Act/Act (ICMA)).
          if ((missing(DCC))|(is.na(DCC))) {
            DCC<-2
            warning("The day count indentifier (DCC) is missing or NA. DCC is set 2 (Act/Act (ICMA))!")
          } else {
            if (!(is.element(DCC,c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)))) {
              DCC<-2
              warning("The day count indentifier (DCC) is not element of {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16}!
                      DCC is set 2 (Act/Act (ICMA))!")
            }
          }
          if ((DCC==7)&((missing(Mat))|(is.na(Mat)))) {
            warning("Maturity date (Mat) is missing or NA. Accrued interest computation for the specified
                    day count convention 30E/360 (ISDA) requires a valid Mat value. NA created!")
          } else {
            if ((DCC==9)&((missing(EOM))|(is.na(EOM)))) {
              warning("End-of-Month-Rule identifier (EOM) is missing or NA. Accrued interest computation for
                      the specified day count convention 30/360 US requires a valid EOM value. NA created!")
            } else {
              if ((DCC==13)&((missing(EOM))|(is.na(EOM)))) {
                warning("End-of-Month-Rule identifier (EOM) is missing or NA. Accrued interest computation for
                      the specified day count convention 30/365 requires a valid EOM value. NA created!")
              } else {
                if ((DCC==2)&((missing(EOM))|(is.na(EOM)))) {
                  warning("End-of-Month-Rule identifier (EOM) is missing or NA. Accrued interest computation for
                          the specified day count convention Act/Act (ICMA) requires a valid EOM value. NA created!")
                } else {
                  if ((DCC==14)&((missing(EOM))|(is.na(EOM)))) {
                    warning("End-of-Month-Rule identifier (EOM) is missing or NA. Accrued interest computation for
                          the specified day count convention Act/365 (Canadian Bond) requires a valid EOM value. NA created!")
                  } else {
                    if ((DCC==4)&(CpY!=1)&((missing(YearNCP))|(is.na(YearNCP)))) {
                      warning(paste0("Year figure of the next coupon payment date after EndDate (YearNCP) is missing or NA.
                                     Accrued interest computation for the specified day count convention Act/365L with ",CpY,
                                     " requires a valid YearNCP value. NA created!"))
                    } else {
                      if ((DCC==2)&((missing(CpY))|(is.na(CpY)))) {
                        CpY<-2
                        warning("Number of interest payments per year (CpY) is missing or NA. Accrued interest computation for
                                the specified day count convention Act/Act (ICMA) requires a valid CpY value. CpY is set 2!")
                      } else {
                        if ((DCC==2)&(!(is.element(CpY,c(1,2,3,4,6,12))))) {
                          CpY<-2
                          warning("Number of interest payments per year (CpY) is not element of {1,2,3,4,6,12}. Accrued interest computation for
                                  the specified day count convention Act/Act (ICMA) requires a valid CpY value. CpY is set 2!")
                        }
                      }
                      if ((DCC==4)&((missing(CpY))|(is.na(CpY)))) {
                        CpY<-2
                        warning("Number of interest payments per year (CpY) is missing or NA. Accrued interest computation for
                                the specified day count convention Act/365L requires a valid CpY value. CpY is set 2!")
                      } else {
                        if ((DCC==4)&(!(is.element(CpY,c(1,2,3,4,6,12))))) {
                          CpY<-2
                          warning("Number of interest payments per year (CpY) is not element of {1,2,3,4,6,12}. Accrued interest computation for
                                  the specified day count convention Act/365L requires a valid CpY value. CpY is set 2!")
                        }
                      }
                      if ((DCC==14)&((missing(CpY))|(is.na(CpY)))) {
                        CpY<-2
                        warning("Number of interest payments per year (CpY) is missing or NA. Accrued interest computation for
                                the specified day count convention Act/365 (Canadian Bond) requires a valid CpY value. CpY is set 2!")
                      } else {
                        if ((DCC==14)&(!(is.element(CpY,c(1,2,3,4,6,12))))) {
                          CpY<-2
                          warning("Number of interest payments per year (CpY) is not element of {1,2,3,4,6,12}. Accrued interest computation for
                                  the specified day count convention Act/365 (Canadian Bond) requires a valid CpY value. CpY is set 2!")
                        }
                      }
                      if (DCC==16) {
                        if ((missing(CpY))|(is.na(CpY))) {
                          CpY<-2
                          warning("Number of interest payments per year (CpY) is missing or NA. Accrued interest computation for
                                the specified day count convention BusDay/252 (Brazilian) requires a valid CpY value. CpY is set 2!")
                        } else {
                          if (!(is.element(CpY,c(1,2,3,4,6,12)))) {
                            CpY<-2
                            warning("Number of interest payments per year (CpY) is not element of {1,2,3,4,6,12}. Accrued interest computation for
                                  the specified day count convention BusDay/252 (Brazilian) requires a valid CpY value. CpY is set 2!")
                          }
                        }
                      }
                      if ((missing(RV))|(is.na(RV))) {
                        RV<-100
                        warning("Redemption value (RV) is missing or NA. RV is set 100!")
                      }
                      if (is.element(DCC,c(1,3,5,6,8,10,11,12,15))) {
                        # c(RV,Coup,DCC,Y1,M1,D1,Y2,M2,D2)
                        Atoms_StartDate<-as.numeric(unlist(strsplit(as.character(StartDate),split = "-")))
                        Atoms_EndDate<-as.numeric(unlist(strsplit(as.character(EndDate),split = "-")))
                        DIST_Output<-DIST(c(DCC,Atoms_StartDate,Atoms_EndDate))
                        NAccr<-DIST_Output[1]
                        AccrInt<-RV*(Coup/100)*DIST_Output[2]
                      } else {
                        if (DCC==16) {
                          # c(RV,Coup,DCC,Y1,M1,D1,Y2,M2,D2,NonBus)
                          NonBus.Start.End<-length(which((NonBusDays.Brazil$Date>=StartDate)&(NonBusDays.Brazil$Date<EndDate)))
                          Atoms_StartDate<-as.numeric(unlist(strsplit(as.character(StartDate),split = "-")))
                          Atoms_EndDate<-as.numeric(unlist(strsplit(as.character(EndDate),split = "-")))
                          DIST_Output<-DIST(c(DCC,Atoms_StartDate,Atoms_EndDate,NonBus.Start.End))
                          NAccr<-DIST_Output[1]
                          AccrInt<-RV*(((1+Coup/100)^(DIST_Output[2]))-1)
                        } else {
                          if (DCC==4) {
                            # c(RV,Coup,DCC,Y1,M1,D1,Y2,M2,D2,YearNCP,CpY)
                            Atoms_StartDate<-as.numeric(unlist(strsplit(as.character(StartDate),split = "-")))
                            Atoms_EndDate<-as.numeric(unlist(strsplit(as.character(EndDate),split = "-")))
                            DIST_Output<-DIST(c(DCC,Atoms_StartDate,Atoms_EndDate,YearNCP,CpY))
                            NAccr<-DIST_Output[1]
                            AccrInt<-RV*(Coup/100)*DIST_Output[2]
                          } else {
                            if (DCC==7) {
                              # c(RV,Coup,DCC,Y1,M1,D1,Y2,M2,D2,YMat,MMat,DMat)
                              Atoms_StartDate<-as.numeric(unlist(strsplit(as.character(StartDate),split = "-")))
                              Atoms_EndDate<-as.numeric(unlist(strsplit(as.character(EndDate),split = "-")))
                              Atoms_Mat<-as.numeric(unlist(strsplit(as.character(Mat),split = "-")))
                              DIST_Output<-DIST(c(DCC,Atoms_StartDate,Atoms_EndDate,Atoms_Mat))
                              NAccr<-DIST_Output[1]
                              AccrInt<-RV*(Coup/100)*DIST_Output[2]
                            } else {
                              if (is.element(DCC,c(9,13))) {
                                #   c(RV,Coup,DCC,Y1,M1,D1,Y2,M2,D2,EOM)
                                Atoms_StartDate<-as.numeric(unlist(strsplit(as.character(StartDate),split = "-")))
                                Atoms_EndDate<-as.numeric(unlist(strsplit(as.character(EndDate),split = "-")))
                                DIST_Output<-DIST(c(DCC,Atoms_StartDate,Atoms_EndDate,EOM))
                                NAccr<-DIST_Output[1]
                                AccrInt<-RV*(Coup/100)*DIST_Output[2]
                              } else {
                                if (DCC==2|DCC==14) {
                                  if (CpY==1) {N_months<-12}
                                  if (CpY==2) {N_months<-6}
                                  if (CpY==3) {N_months<-4}
                                  if (CpY==4) {N_months<-3}
                                  if (CpY==6) {N_months<-2}
                                  if (CpY==12) {N_months<-1}
                                  AtomVector_StartDate<-as.numeric(unlist(strsplit(as.character(StartDate),split = "-")))
                                  Atom1StartDate<-AtomVector_StartDate[1]
                                  Atom2StartDate<-AtomVector_StartDate[2]
                                  Atom3StartDate<-AtomVector_StartDate[3]
                                  StartDateHelp<-as.Date(paste(Atom1StartDate,Atom2StartDate,15,sep="-"))
                                  AtomVector_EndDate<-as.numeric(unlist(strsplit(as.character(EndDate),split = "-")))
                                  Atom1EndDate<-AtomVector_EndDate[1]
                                  Atom2EndDate<-AtomVector_EndDate[2]
                                  Atom3EndDate<-AtomVector_EndDate[3]
                                  LDM_EndDate<-as.numeric(Date_LDM(c(Atom1EndDate,Atom2EndDate,Atom3EndDate)))
                                  LDM_EndDate<-as.Date(paste(LDM_EndDate[1],LDM_EndDate[2],LDM_EndDate[3],sep="-"))
                                  AnnivDates<-seq(StartDateHelp,LDM_EndDate,by=paste(N_months," months",sep=""))
                                  # assigning the reference date that determines the day figures of all AnnivDates
                                  Atom1Refer<-Atom1StartDate
                                  Atom2Refer<-Atom2StartDate
                                  Atom3Refer<-Atom3StartDate
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
                                  # creating the anniversary date preceding AD1
                                  AtomVector_AD1<-as.numeric(unlist(strsplit(as.character(AnnivDates[1]),split = "-")))
                                  Atom1AD1<-AtomVector_AD1[1]
                                  Atom2AD1<-AtomVector_AD1[2]
                                  Atom3AD1<-AtomVector_AD1[3]
                                  PrevDate<-as.numeric(CppPrevDate(c(Atom1AD1,Atom2AD1,Atom3AD1,Atom1AD1,Atom2AD1,Atom3AD1,Atom1Refer,Atom2Refer,Atom3Refer,CpY,EOM)))
                                  PrevDate<-as.Date(paste(PrevDate[1],PrevDate[2],PrevDate[3],sep="-"))
                                  AnnivDates<-c(PrevDate,AnnivDates)
                                  AnnivDates<-sort(na.omit(AnnivDates[!duplicated(AnnivDates)]))
                                  # creating the anniversary date succeeding ADfin
                                  AtomVector_ADfin<-as.numeric(unlist(strsplit(as.character(AnnivDates[length(AnnivDates)]),split = "-")))
                                  Atom1ADfin<-AtomVector_ADfin[1]
                                  Atom2ADfin<-AtomVector_ADfin[2]
                                  Atom3ADfin<-AtomVector_ADfin[3]
                                  SuccDate<-as.numeric(CppSuccDate(c(Atom1ADfin,Atom2ADfin,Atom3ADfin,Atom1ADfin,Atom2ADfin,Atom3ADfin,Atom1Refer,Atom2Refer,Atom3Refer,CpY,EOM)))
                                  SuccDate<-as.Date(paste(SuccDate[1],SuccDate[2],SuccDate[3],sep="-"))
                                  AnnivDates<-c(AnnivDates,SuccDate)
                                  AnnivDates<-sort(na.omit(AnnivDates[!duplicated(AnnivDates)]))
                                  # DIST: for DCC = 2 x is a vector of 22 integers:
                                  #   c(DCC,Y1,M1,D1,Y2,M2,D2,Y3,M3,D3,Y4,M4,D4,Y5,M5,D5,Y6,M6,D6,P,N,CpY) with
                                  # Y1-M1-D1 = PCD(t_a,AD) ; Y2-M2-D2 = t_a ; Y3-M3-D3 = NCD(t_a,AD)
                                  # Y4-M4-D4 = PCD(t_b,AD) ; Y5-M5-D5 = t_b ; Y6-M6-D6 = NCD(t_b,AD)
                                  # P = P(t_b,AD) ; N = N(t_a,AD)
                                  PCD_StartDate<-PCD(StartDate,AnnivDates)
                                  Atoms_PCD_StartDate<-as.numeric(unlist(strsplit(as.character(PCD_StartDate),split = "-")))
                                  # Atom1_PCD_StartDate<-Atoms_PCD_StartDate[1]
                                  # Atom2_PCD_StartDate<-Atoms_PCD_StartDate[2]
                                  # Atom3_PCD_StartDate<-Atoms_PCD_StartDate[3]
                                  NCD_StartDate<-NCD(StartDate,AnnivDates)
                                  Atoms_NCD_StartDate<-as.numeric(unlist(strsplit(as.character(NCD_StartDate),split = "-")))
                                  # Atom1_NCD_StartDate<-Atoms_NCD_StartDate[1]
                                  # Atom2_NCD_StartDate<-Atoms_NCD_StartDate[2]
                                  # Atom3_NCD_StartDate<-Atoms_NCD_StartDate[3]
                                  PCD_EndDate<-PCD(EndDate,AnnivDates)
                                  Atoms_PCD_EndDate<-as.numeric(unlist(strsplit(as.character(PCD_EndDate),split = "-")))
                                  # Atom1_PCD_EndDate<-Atoms_PCD_EndDate[1]
                                  # Atom2_PCD_EndDate<-Atoms_PCD_EndDate[2]
                                  # Atom3_PCD_EndDate<-Atoms_PCD_EndDate[3]
                                  NCD_EndDate<-NCD(EndDate,AnnivDates)
                                  Atoms_NCD_EndDate<-as.numeric(unlist(strsplit(as.character(NCD_EndDate),split = "-")))
                                  # Atom1_NCD_EndDate<-Atoms_NCD_EndDate[1]
                                  # Atom2_NCD_EndDate<-Atoms_NCD_EndDate[2]
                                  # Atom3_NCD_EndDate<-Atoms_NCD_EndDate[3]
                                  AD_indexes<-c(1:length(AnnivDates))-1
                                  AD_List<-list(AnnivDates,AD_indexes)
                                  N_StartDate<-AD_List[[2]][which(AD_List[[1]]==NCD_StartDate)]
                                  P_EndDate<-AD_List[[2]][which(AD_List[[1]]==PCD_EndDate)]
                                  DIST_Output<-DIST(c(DCC,Atoms_PCD_StartDate,AtomVector_StartDate,Atoms_NCD_StartDate,
                                                      Atoms_PCD_EndDate,AtomVector_EndDate,Atoms_NCD_EndDate,P_EndDate,N_StartDate,CpY))
                                  NAccr<-DIST_Output[1]
                                  AccrInt<-RV*(Coup/100)*DIST_Output[2]
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  Output<-list(AccrInt=AccrInt,DaysAccrued=NAccr)
  return(Output)
}




