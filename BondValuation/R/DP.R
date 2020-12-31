#' DP (dirty price calculation of a fixed-coupon bond)
#'
#' \bold{DP} returns a bond's temporal and pecuniary characteristics on the desired calendar date
#' according to the methodology presented in Djatschenko (2018).
#'
#' The function \bold{DP} generates a list of the two data frames \code{Dates} and \code{Cash},
#' which contain the relevant date-related and pecuniary characteristics that were either provided
#' by the user or calculated by the function. \bold{Value} provides further information on the
#' output.
#'
#' @param CP The bond's clean price.
#' @param SETT The settlement date. Date class object with format "\%Y-\%m-\%d". (required)
#' @param Em The bond's issue date. Date class object with format "\%Y-\%m-\%d". (required)
#' @param Mat So-called "maturity date" i.e. date on which the redemption value and the final interest
#'        are paid. Date class object with format "\%Y-\%m-\%d". (required)
#' @param CpY Number of interest payments per year (non-negative integer; element of the set
#'        \{0,1,2,3,4,6,12\}. Default: 2.
#' @param FIPD First interest payment date after \code{Em}. Date class object with format "\%Y-\%m-\%d". Default: \code{NA}.
#' @param LIPD Last interest payment date before \code{Mat}. Date class object with format "\%Y-\%m-\%d". Default: \code{NA}.
#' @param FIAD Date on which the interest accrual starts (so-called "dated date"). Date class object with format "\%Y-\%m-\%d". Default: \code{NA}.
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
#' @param AnnivDatesOutput A list containing the output of the function AnnivDates. Default: \code{NA}.
#'
#'
#' @return
#'   \describe{
#'     \item{\emph{\bold{Dates}} (data frame)}{
#'       \describe{
#'         \item{\emph{Previous_CouponDate}}{}
#'         \item{\emph{SettlementDate}}{}
#'         \item{\emph{Next_CouponDate}}{}
#'         \item{\emph{DaysAccrued}}{The number of days accrued from \emph{Previous_CouponDate} to
#'                                   \emph{Next_CouponDate}, incl. the earlier and excl. the later
#'                                   date.}
#'         \item{\emph{DaysInPeriod}}{The number of interest accruing days in the coupon period
#'                                    from \emph{Previous_CouponDate} to \emph{Next_CouponDate}.}
#'       }
#'     }
#'     \item{\emph{\bold{Cash}} (data frame)}{
#'       \describe{
#'         \item{\emph{Dirty_Price}}{Sum of \emph{Clean_Price} and \emph{Accrued_Interest}.}
#'         \item{\emph{Clean_Price}}{The clean price entered.}
#'         \item{\emph{Accrued_Interest}}{The amount of accrued interest on \emph{SettlementDate}.}
#'         \item{\emph{CouponPayment}}{The interest payment on \emph{Next_CouponDate}.}
#'       }
#'     }
#'   }
#'
#' @references
#' \enumerate{
#'   \item{Djatschenko, Wadim, The Nitty Gritty of Bond Valuation: A Generalized Methodology for Fixed Coupon Bond Analysis Allowing for Irregular Periods and Various Day Count Conventions (November 5, 2018). Available at SSRN: https://ssrn.com/abstract=3205167.}
#' }
#'
#' @examples
#' CP<-rep(100,16)
#' SETT<-rep(as.Date("2014-10-15"),16)
#' Em<-rep(as.Date("2013-11-30"),16)
#' Mat<-rep(as.Date("2021-04-21"),16)
#' CpY<-rep(2,16)
#' FIPD<-rep(as.Date("2015-02-28"),16)
#' LIPD<-rep(as.Date("2020-02-29"),16)
#' FIAD<-rep(as.Date("2013-11-30"),16)
#' RV<-rep(100,16)
#' Coup<-rep(5.25,16)
#' DCC<-seq(1,16,by=1)
#' DP.DCC_Comparison<-data.frame(CP,SETT,Em,Mat,CpY,FIPD,LIPD,FIAD,RV,Coup,DCC)
#'
#' # you can pass an array to AnnivDates
#' List<-suppressWarnings(
#'         AnnivDates(unlist(DP.DCC_Comparison[1,c(3:11)],use.names=FALSE))
#' )
#'
#' # and use its output in DP
#' suppressWarnings(
#'        DP(unlist(DP.DCC_Comparison[1,c(1:11)],use.names=FALSE),AnnivDatesOutput=List)
#' )
#'
#' # or just apply DP to the data frame
#' DP.Output<-suppressWarnings(
#'               apply(DP.DCC_Comparison[,c('CP','SETT','Em','Mat','CpY','FIPD',
#'                                            'LIPD','FIAD','RV','Coup','DCC')],
#'                      1,function(y) DP(y[1],y[2],y[3],y[4],y[5],y[6],y[7],
#'                                       y[8],y[9],y[10],y[11])))
#'
#' DiryPrice<-do.call(rbind,lapply(lapply(DP.Output, `[[`, 2), `[[`, 1))
#' DP.DCC_Comparison<-cbind(DP.DCC_Comparison,DiryPrice)
#' DP.DCC_Comparison
#'
#' @export
#'
DP<-function(CP=as.numeric(NA),SETT=as.Date(NA),Em=as.Date(NA),Mat=as.Date(NA),CpY=as.numeric(NA),FIPD=as.Date(NA),LIPD=as.Date(NA),FIAD=as.Date(NA),RV=as.numeric(NA),Coup=as.numeric(NA),DCC=as.numeric(NA),EOM=as.numeric(NA),DateOrigin=as.Date("1970-01-01"),InputCheck=1,FindEOM=FALSE,RegCF.equal=0,AnnivDatesOutput=as.list(NA)) {
  if (length(CP)>1) {
    arglist<-CP
    argnames<-c("CP","SETT","Em","Mat","CpY","FIPD","LIPD","FIAD","RV","Coup","DCC","EOM","DateOrigin","InputCheck","FindEOM","RegCF.equal")
    for (i in c(1:length(arglist))) {
      assign(argnames[i],arglist[i])
    }
  }
  if (InputCheck==1) {
    CheckedInput<-InputFormatCheck(CP=CP,SETT=SETT,Em=Em,Mat=Mat,CpY=CpY,FIPD=FIPD,LIPD=LIPD,FIAD=FIAD,RV=RV,Coup=Coup,DCC=DCC,EOM=EOM,DateOrigin=DateOrigin)
    CP<-CheckedInput$CP
    SETT<-CheckedInput$SETT
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
  DP<-as.numeric(NA)
  if ((missing(CP))|(is.na(CP))) {
    CP<-0
    warning("Clean price (CP) is missing or NA. CP is set 0.")
  }
  AccrInt<-as.numeric(NA)
  CouponPayment<-as.numeric(NA)
  NAccr<-as.numeric(NA)
  NPeriod<-as.numeric(NA)
  PCD_SETT<-as.Date(NA)
  NCD_SETT<-as.Date(NA)
  if ((missing(SETT))|(is.na(SETT))) {
    SETT<-as.Date(NA)
    warning("Settlement date (SETT) is missing or NA. NA created!")
  } else {
    if (all(is.na(AnnivDatesOutput))) {
      BondAnalysis<-suppressWarnings(AnnivDates(Em=Em,Mat=Mat,CpY=CpY,FIPD=FIPD,LIPD=LIPD,FIAD=FIAD,RV=RV,Coup=Coup,DCC=DCC,EOM=EOM,DateOrigin=DateOrigin,InputCheck=0,FindEOM=FindEOM,RegCF.equal=RegCF.equal))
    } else {
      BondAnalysis<-AnnivDatesOutput
    }
    RealDates<-na.omit(BondAnalysis$DateVectors$RealDates)
    AnnivDates<-na.omit(BondAnalysis$DateVectors$AnnivDates)
    if (is.na(BondAnalysis$Traits$FIPD)) {
      FIPD<-BondAnalysis$Traits$est_FIPD
    } else {
      FIPD<-BondAnalysis$Traits$FIPD
    }
    if (is.na(BondAnalysis$Traits$LIPD)) {
      LIPD<-BondAnalysis$Traits$est_LIPD
    } else {
      LIPD<-BondAnalysis$Traits$LIPD
    }
    if (BondAnalysis$Warnings$EmMatMissing==1) {
      warning("Maturity date (Mat) is missing or NA. NA created!")
    } else {
      if (BondAnalysis$Warnings$NegLifeFlag==1) {
        warning("Issue date (Em) is not before maturity date (Mat)! NA created!")
      } else {
        if ((missing(Coup))|(is.na(Coup))) {
          warning("The supplied interest rate p.a. (Coup) is NA or cannot be processed. NA created!")
        } else {
          if ((SETT<RealDates[1])|(SETT>=RealDates[length(RealDates)])) {
            warning("Settlement date (SETT) is not between issue date (Em) and maturity date (Mat). NA created!")
          } else {
            if (BondAnalysis$Warnings$ZeroFlag==1) {
              DP<-CP
              AccrInt<-as.numeric(0)
              NAccr<-as.numeric(NA)
              NPeriod<-as.numeric(NA)
              PCD_SETT<-Em
              NCD_SETT<-Mat
              CouponPayment<-BondAnalysis$PaySched$CoupPayments[1]
              warning("This is a Zero Coupon bond! No interest accrues!")
            } else {
              EOM<-BondAnalysis$Traits$EOM_used
              # If DCC is not provided or NA or not element of {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16}, the following code sets it 2 (Act/Act (ICMA)).
              if ((missing(DCC))|(is.na(DCC))) {
                DCC<-2
                warning("The day count indentifier (DCC) is missing or NA. DCC is set 2 (Act/Act (ICMA))!")
              } else {
                if (!(is.element(DCC,c(1:16)))) {
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
                          PCD_SETT<-PCD(SETT,RealDates)
                          Atoms_PCD_SETT<-as.numeric(unlist(strsplit(as.character(PCD_SETT),split = "-")))
                          Atoms_SETT<-as.numeric(unlist(strsplit(as.character(SETT),split = "-")))
                          NCD_SETT<-NCD(SETT,RealDates)
                          Atoms_NCD_SETT<-as.numeric(unlist(strsplit(as.character(NCD_SETT),split = "-")))
                          DIST_Output_AccrInt<-DIST(c(DCC,Atoms_PCD_SETT,Atoms_SETT))
                          DIST_Output_Coup<-DIST(c(DCC,Atoms_PCD_SETT,Atoms_NCD_SETT))
                          NAccr<-DIST_Output_AccrInt[1]
                          AccrInt<-RV*(Coup/100)*DIST_Output_AccrInt[2]
                          NPeriod<-DIST_Output_Coup[1]
                          CouponPayment<-RV*(Coup/100)*DIST_Output_Coup[2]
                          DP<-CP+AccrInt
                        } else {
                          if (DCC==16) {
                            # c(RV,Coup,DCC,Y1,M1,D1,Y2,M2,D2,NonBus)
                            PCD_SETT<-PCD(SETT,RealDates)
                            Atoms_PCD_SETT<-as.numeric(unlist(strsplit(as.character(PCD_SETT),split = "-")))
                            Atoms_SETT<-as.numeric(unlist(strsplit(as.character(SETT),split = "-")))
                            NCD_SETT<-NCD(SETT,RealDates)
                            Atoms_NCD_SETT<-as.numeric(unlist(strsplit(as.character(NCD_SETT),split = "-")))
                            NonBus.PCD.SETT<-length(which((NonBusDays.Brazil$Date>=PCD_SETT)&(NonBusDays.Brazil$Date<SETT)))
                            NonBus.PCD.NCD<-length(which((NonBusDays.Brazil$Date>=PCD_SETT)&(NonBusDays.Brazil$Date<NCD_SETT)))
                            DIST_Output_AccrInt<-DIST(c(DCC,Atoms_PCD_SETT,Atoms_SETT,NonBus.PCD.SETT))
                            DIST_Output_Coup<-DIST(c(DCC,Atoms_PCD_SETT,Atoms_NCD_SETT,NonBus.PCD.NCD))
                            if ((!(SETT<FIPD))&(!(SETT>LIPD))) {
                              CouponPayment<-RV*(((1+Coup/100)^(1/CpY))-1)
                            } else {
                              CouponPayment<-RV*(((1+Coup/100)^(DIST_Output_Coup[2]))-1)
                            }
                            NAccr<-DIST_Output_AccrInt[1]
                            AccrInt<-RV*(((1+Coup/100)^(DIST_Output_AccrInt[2]))-1)
                            if (AccrInt>CouponPayment) {
                              AccrInt<-CouponPayment
                            }
                            NPeriod<-DIST_Output_Coup[1]
                            DP<-CP+AccrInt
                          } else {
                            if (DCC==4) {
                              # c(RV,Coup,DCC,Y1,M1,D1,Y2,M2,D2,YearNCP,CpY)
                              PCD_SETT<-PCD(SETT,RealDates)
                              Atoms_PCD_SETT<-as.numeric(unlist(strsplit(as.character(PCD_SETT),split = "-")))
                              Atoms_SETT<-as.numeric(unlist(strsplit(as.character(SETT),split = "-")))
                              NCD_SETT<-NCD(SETT,RealDates)
                              Atoms_NCD_SETT<-as.numeric(unlist(strsplit(as.character(NCD_SETT),split = "-")))
                              DIST_Output_AccrInt<-DIST(c(DCC,Atoms_PCD_SETT,Atoms_SETT,Atoms_NCD_SETT[1],CpY))
                              DIST_Output_Coup<-DIST(c(DCC,Atoms_PCD_SETT,Atoms_NCD_SETT,Atoms_NCD_SETT[1],CpY))
                              NAccr<-DIST_Output_AccrInt[1]
                              AccrInt<-RV*(Coup/100)*DIST_Output_AccrInt[2]
                              NPeriod<-DIST_Output_Coup[1]
                              CouponPayment<-RV*(Coup/100)*DIST_Output_Coup[2]
                              DP<-CP+AccrInt
                            } else {
                              if (DCC==7) {
                                # c(RV,Coup,DCC,Y1,M1,D1,Y2,M2,D2,YMat,MMat,DMat)
                                PCD_SETT<-PCD(SETT,RealDates)
                                Atoms_PCD_SETT<-as.numeric(unlist(strsplit(as.character(PCD_SETT),split = "-")))
                                Atoms_SETT<-as.numeric(unlist(strsplit(as.character(SETT),split = "-")))
                                NCD_SETT<-NCD(SETT,RealDates)
                                Atoms_NCD_SETT<-as.numeric(unlist(strsplit(as.character(NCD_SETT),split = "-")))
                                Atoms_Mat<-as.numeric(unlist(strsplit(as.character(Mat),split = "-")))
                                DIST_Output_AccrInt<-DIST(c(DCC,Atoms_PCD_SETT,Atoms_SETT,Atoms_Mat))
                                DIST_Output_Coup<-DIST(c(DCC,Atoms_PCD_SETT,Atoms_NCD_SETT,Atoms_Mat))
                                NAccr<-DIST_Output_AccrInt[1]
                                AccrInt<-RV*(Coup/100)*DIST_Output_AccrInt[2]
                                NPeriod<-DIST_Output_Coup[1]
                                CouponPayment<-RV*(Coup/100)*DIST_Output_Coup[2]
                                DP<-CP+AccrInt
                              } else {
                                if (is.element(DCC,c(9,13))) {
                                  #   c(RV,Coup,DCC,Y1,M1,D1,Y2,M2,D2,EOM)
                                  PCD_SETT<-PCD(SETT,RealDates)
                                  Atoms_PCD_SETT<-as.numeric(unlist(strsplit(as.character(PCD_SETT),split = "-")))
                                  Atoms_SETT<-as.numeric(unlist(strsplit(as.character(SETT),split = "-")))
                                  NCD_SETT<-NCD(SETT,RealDates)
                                  Atoms_NCD_SETT<-as.numeric(unlist(strsplit(as.character(NCD_SETT),split = "-")))
                                  DIST_Output_AccrInt<-DIST(c(DCC,Atoms_PCD_SETT,Atoms_SETT,EOM))
                                  DIST_Output_Coup<-DIST(c(DCC,Atoms_PCD_SETT,Atoms_NCD_SETT,EOM))
                                  NAccr<-DIST_Output_AccrInt[1]
                                  AccrInt<-RV*(Coup/100)*DIST_Output_AccrInt[2]
                                  NPeriod<-DIST_Output_Coup[1]
                                  CouponPayment<-RV*(Coup/100)*DIST_Output_Coup[2]
                                  DP<-CP+AccrInt
                                } else {
                                  if (DCC==2|DCC==14) {
                                    PCD_SETT<-PCD(SETT,RealDates)
                                    NCD_SETT<-NCD(SETT,RealDates)
                                    Refer<-BondAnalysis$Traits$Refer
                                    AD_indexes<-BondAnalysis$DateVectors$AD_indexes
                                    AD_indexes<-c((AD_indexes[1]-1),AD_indexes,(AD_indexes[length(AD_indexes)]+1))
                                    AtomVector_Refer<-as.numeric(unlist(strsplit(as.character(Refer),split = "-")))
                                    Atom1Refer<-AtomVector_Refer[1]
                                    Atom2Refer<-AtomVector_Refer[2]
                                    Atom3Refer<-AtomVector_Refer[3]
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
                                    AD_List<-list(AnnivDates,AD_indexes)
                                    # DIST: for DCC = {2,14} x is a vector of 22 integers:
                                    #   c(DCC,Y1,M1,D1,Y2,M2,D2,Y3,M3,D3,Y4,M4,D4,Y5,M5,D5,Y6,M6,D6,P,N,CpY) with
                                    # Y1-M1-D1 = PCD(t_a,AD) ; Y2-M2-D2 = t_a ; Y3-M3-D3 = NCD(t_a,AD)
                                    # Y4-M4-D4 = PCD(t_b,AD) ; Y5-M5-D5 = t_b ; Y6-M6-D6 = NCD(t_b,AD)
                                    # P = P(t_b,AD) ; N = N(t_a,AD)
                                    # for AccrInt: t_a = PCD(SETT,SD) and t_b = SETT
                                    ta_AccrInt<-PCD(SETT,RealDates)
                                    Atoms_ta_AccrInt<-as.numeric(unlist(strsplit(as.character(ta_AccrInt),split = "-")))
                                    tb_AccrInt<-SETT
                                    Atoms_tb_AccrInt<-as.numeric(unlist(strsplit(as.character(tb_AccrInt),split = "-")))
                                    PCD_ta_AccrInt<-PCD(ta_AccrInt,AnnivDates)
                                    Atoms_PCD_ta_AccrInt<-as.numeric(unlist(strsplit(as.character(PCD_ta_AccrInt),split = "-")))
                                    NCD_ta_AccrInt<-NCD(ta_AccrInt,AnnivDates)
                                    Atoms_NCD_ta_AccrInt<-as.numeric(unlist(strsplit(as.character(NCD_ta_AccrInt),split = "-")))
                                    PCD_tb_AccrInt<-PCD(tb_AccrInt,AnnivDates)
                                    Atoms_PCD_tb_AccrInt<-as.numeric(unlist(strsplit(as.character(PCD_tb_AccrInt),split = "-")))
                                    NCD_tb_AccrInt<-NCD(tb_AccrInt,AnnivDates)
                                    Atoms_NCD_tb_AccrInt<-as.numeric(unlist(strsplit(as.character(NCD_tb_AccrInt),split = "-")))
                                    N_ta_AccrInt<-AD_List[[2]][which(AD_List[[1]]==NCD_ta_AccrInt)]
                                    P_tb_AccrInt<-AD_List[[2]][which(AD_List[[1]]==PCD_tb_AccrInt)]
                                    DIST_Output_AccrInt<-DIST(c(DCC,Atoms_PCD_ta_AccrInt,Atoms_ta_AccrInt,Atoms_NCD_ta_AccrInt,
                                                                Atoms_PCD_tb_AccrInt,Atoms_tb_AccrInt,Atoms_NCD_tb_AccrInt,P_tb_AccrInt,N_ta_AccrInt,CpY))
                                    NAccr<-DIST_Output_AccrInt[1]
                                    AccrInt<-RV*(Coup/100)*DIST_Output_AccrInt[2]
                                    # for Coup: t_a = PCD(SETT,SD) and t_b = NCD(SETT,SD)
                                    ta_Coup<-PCD(SETT,RealDates)
                                    Atoms_ta_Coup<-as.numeric(unlist(strsplit(as.character(ta_Coup),split = "-")))
                                    tb_Coup<-NCD(SETT,RealDates)
                                    Atoms_tb_Coup<-as.numeric(unlist(strsplit(as.character(tb_Coup),split = "-")))
                                    PCD_ta_Coup<-PCD(ta_Coup,AnnivDates)
                                    Atoms_PCD_ta_Coup<-as.numeric(unlist(strsplit(as.character(PCD_ta_Coup),split = "-")))
                                    NCD_ta_Coup<-NCD(ta_Coup,AnnivDates)
                                    Atoms_NCD_ta_Coup<-as.numeric(unlist(strsplit(as.character(NCD_ta_Coup),split = "-")))
                                    PCD_tb_Coup<-PCD(tb_Coup,AnnivDates)
                                    Atoms_PCD_tb_Coup<-as.numeric(unlist(strsplit(as.character(PCD_tb_Coup),split = "-")))
                                    NCD_tb_Coup<-NCD(tb_Coup,AnnivDates)
                                    Atoms_NCD_tb_Coup<-as.numeric(unlist(strsplit(as.character(NCD_tb_Coup),split = "-")))
                                    N_ta_Coup<-AD_List[[2]][which(AD_List[[1]]==NCD_ta_Coup)]
                                    P_tb_Coup<-AD_List[[2]][which(AD_List[[1]]==PCD_tb_Coup)]
                                    DIST_Output_Coup<-DIST(c(DCC,Atoms_PCD_ta_Coup,Atoms_ta_Coup,Atoms_NCD_ta_Coup,
                                                             Atoms_PCD_tb_Coup,Atoms_tb_Coup,Atoms_NCD_tb_Coup,P_tb_Coup,N_ta_Coup,CpY))
                                    NPeriod<-DIST_Output_Coup[1]
                                    CouponPayment<-RV*(Coup/100)*DIST_Output_Coup[2]
                                    DP<-CP+AccrInt
                                  }
                                }
                              }
                            }
                          }
                        }
                        if ((!(RegCF.equal==0))&(!(DCC==16))) {
                          if ((!(SETT<FIPD))&(!(SETT>LIPD))) {
                            CouponPayment<-RV*(Coup/(CpY*100))
                            if (AccrInt>CouponPayment) {
                              AccrInt<-CouponPayment
                              DP<-CP+AccrInt
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
  Dates<-data.frame(Previous_CouponDate=PCD_SETT,SettlementDate=SETT,Next_CouponDate=NCD_SETT,DaysAccrued=NAccr,DaysInPeriod=NPeriod)
  Cash<-data.frame(Dirty_Price=DP,Clean_Price=CP,Accrued_Interest=AccrInt,CouponPayment=CouponPayment)
  DP_Out<-list(Dates,Cash)
  return(DP_Out)
}



