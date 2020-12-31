#' BondVal.Price (calculation of CP, AccrInt, DP, ModDUR, MacDUR and Conv)
#'
#' \bold{BondVal.Price} computes a bond's clean price given its yield.
#'
#' The function \bold{BondVal.Price} uses the function \bold{AnnivDates} to analyze the bond
#' and computes the clean price, the accrued interest, the dirty price and the sensitivity
#' measures modified duration (ModDUR), MacAulay duration (MacDUR) and convexity according
#' to the methodology presented in Djatschenko (2018).
#'
#' @param YtM The bond's yield to maturity p.a. on \code{SETT}. (required)
#' @param SETT The settlement date. Date class object with format "\%Y-\%m-\%d". (required)
#' @param Em The bond's issue date. Date class object with format "\%Y-\%m-\%d". (required)
#' @param Mat So-called "maturity date" i.e. date on which the redemption value and the final interest
#'        are paid. Date class object with format "\%Y-\%m-\%d". (required)
#' @param CpY Number of interest payments per year (non-negative integer; element of the set
#'        \{0,1,2,3,4,6,12\}. Default: 2.
#' @param FIPD First interest payment date after \code{Em}. Date class object with format "\%Y-\%m-\%d". Default: \code{NA}.
#' @param LIPD Last interest payment date before \code{Mat}. Date class object with format "\%Y-\%m-\%d". Default: \code{NA}.
#' @param FIAD Date on which the interest accrual starts (so-called "dated date"). Date class object with format "\%Y-\%m-\%d". Default: \code{NA}.
#' @param RV The redemption value of the bond. Default: \code{100}.
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
#' @param SimpleLastPeriod Specifies the interest calculation method in the final coupon period. Default: \code{TRUE}.
#' @param Calc.Method If 1, discount powers are computed with the same DCC as accrued interest.
#'        If 0, discount powers are computed with DCC=2. Default: 1.
#' @param AnnivDatesOutput A list containing the output of the function AnnivDates. Default: \code{NA}.
#'
#' @return
#'   \describe{
#'     \item{CP}{The bond's clean price.}
#'     \item{AccrInt}{The amount of accrued interest.}
#'     \item{DP}{The bond's dirty price.}
#'     \item{ytm.p.a.}{Annualized yield to maturity.}
#'     \item{ModDUR.inYears}{Modified duration in years.}
#'     \item{MacDUR.inYears}{MacAulay duration in years.}
#'     \item{Conv.inYears}{Convexity in years.}
#'     \item{ModDUR.inPeriods}{Modified duration in periods.}
#'     \item{MacDUR.inPeriods}{MacAulay duration in periods.}
#'     \item{Conv.inPeriods}{Convexity in periods.}
#'     \item{tau}{Relative Position of the settlement date in regular periods.}
#'   }
#'
#' @references
#' \enumerate{
#'   \item{Djatschenko, Wadim, The Nitty Gritty of Bond Valuation: A Generalized Methodology for Fixed Coupon Bond Analysis Allowing for Irregular Periods and Various Day Count Conventions (November 5, 2018). Available at SSRN: https://ssrn.com/abstract=3205167.}
#' }
#'
#' @examples
#' data(PanelSomeBonds2016)
#' randombond<-sample(c(1:length(which(!(duplicated(PanelSomeBonds2016$ID.No))))),1)
#' df.randombond<-PanelSomeBonds2016[which(PanelSomeBonds2016$ID.No==randombond),]
#'
#' PreAnalysis.randombond<-suppressWarnings(AnnivDates(
#'   unlist(df.randombond[
#'            1,c('Issue.Date','Mat.Date','CpY.Input','FIPD.Input','LIPD.Input',
#'                'FIAD.Input','RV.Input','Coup.Input','DCC.Input','EOM.Input')],
#'          use.names=FALSE)))
#'
#' system.time(
#'   for (i in c(1:nrow(df.randombond))) {
#'     BondVal.Price.Output<-suppressWarnings(BondVal.Price(
#'       unlist(
#'         df.randombond[
#'           i,c('YtM.Input','TradeDate','Issue.Date','Mat.Date','CpY.Input',
#'               'FIPD.Input','LIPD.Input','FIAD.Input','RV.Input','Coup.Input',
#'               'DCC.Input','EOM.Input')],use.names=FALSE),
#'       AnnivDatesOutput=PreAnalysis.randombond))
#'     df.randombond$CP.Out[i]<-BondVal.Price.Output$CP
#'   }
#' )
#' plot(seq(1,nrow(df.randombond),by=1),df.randombond$CP.Out,"l")
#'
#' @export
#'
BondVal.Price<-function(YtM=as.numeric(NA),SETT=as.Date(NA),Em=as.Date(NA),Mat=as.Date(NA),CpY=as.numeric(NA),FIPD=as.Date(NA),LIPD=as.Date(NA),FIAD=as.Date(NA),RV=as.numeric(NA),Coup=as.numeric(NA),DCC=as.numeric(NA),EOM=as.numeric(NA),DateOrigin=as.Date("1970-01-01"),InputCheck=1,FindEOM=FALSE,RegCF.equal=0,SimpleLastPeriod=TRUE,Calc.Method=1,AnnivDatesOutput=as.list(NA)) {
  if (length(YtM)>1) {
    arglist<-YtM
    argnames<-c("YtM","SETT","Em","Mat","CpY","FIPD","LIPD","FIAD","RV","Coup","DCC","EOM","DateOrigin","InputCheck","FindEOM","RegCF.equal","SimpleLastPeriod","Calc.Method")
    for (i in c(1:length(arglist))) {
      assign(argnames[i],arglist[i])
    }
  }
  if (InputCheck==1) {
    CheckedInput<-InputFormatCheck(YtM=YtM,SETT=SETT,Em=Em,Mat=Mat,CpY=CpY,FIPD=FIPD,LIPD=LIPD,FIAD=FIAD,RV=RV,Coup=Coup,DCC=DCC,EOM=EOM,DateOrigin=DateOrigin)
    YtM<-CheckedInput$YtM
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
  CP<-as.numeric(NA)
  n<-as.numeric(NA)
  k<-as.numeric(NA)
  CN_tau<-as.numeric(NA)
  tau<-as.numeric(NA)
  w<-as.numeric(NA)
  eta<-as.numeric(NA)
  z<-as.numeric(NA)
  AccrInt<-as.numeric(NA)
  ModDUR.inYears<-as.numeric(NA)
  MacDUR.inYears<-as.numeric(NA)
  Conv.inYears<-as.numeric(NA)
  ModDUR.inPeriods<-as.numeric(NA)
  MacDUR.inPeriods<-as.numeric(NA)
  Conv.inPeriods<-as.numeric(NA)
  CF_final<-as.numeric(NA)
  CF.remain<-as.numeric(NA)
  DiscPowerVector<-as.numeric(NA)
  PriceEqn<-as.numeric(NA)
  if ((missing(YtM))|(is.na(YtM))) {
    YtM<-as.numeric(NA)
    warning("Yield to maturity (YtM) is missing or NA. NA created!")
  } else {
    if ((missing(SETT))|(is.na(SETT))) {
      SETT<-as.Date(NA)
      warning("Settlement date (SETT) is missing or NA. NA created!")
    } else {
      if ((SETT<Em)|(Mat<=SETT)) {
        warning("Settlement date (SETT) is outside bond's lifespan. NA created!")
      } else {
        if (Calc.Method==0) {
          RegCF.equal<-1
          BondAnalysis<-suppressWarnings(AnnivDates(Em=Em,Mat=Mat,CpY=CpY,FIPD=FIPD,LIPD=LIPD,FIAD=FIAD,RV=RV,Coup=Coup,DCC=DCC,EOM=EOM,DateOrigin=DateOrigin,InputCheck=0,FindEOM=FindEOM,RegCF.equal=RegCF.equal))
        } else {
          if (all(is.na(AnnivDatesOutput))) {
            BondAnalysis<-suppressWarnings(AnnivDates(Em=Em,Mat=Mat,CpY=CpY,FIPD=FIPD,LIPD=LIPD,FIAD=FIAD,RV=RV,Coup=Coup,DCC=DCC,EOM=EOM,DateOrigin=DateOrigin,InputCheck=0,FindEOM=FindEOM,RegCF.equal=RegCF.equal))
          } else {
            BondAnalysis<-AnnivDatesOutput
          }
        }
        CpY<-BondAnalysis$Traits$CpY
        RV<-BondAnalysis$Traits$Par
        RealDates<-na.omit(BondAnalysis$DateVectors$RealDates)
        if ((SETT<RealDates[1])|(SETT>=RealDates[length(RealDates)])) {
          warning("Settlement date (SETT) is not between issue date (Em) and maturity date (Mat). NA created!")
        } else {
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
          DP.Output<-DP(RV,SETT,Em,Mat,CpY,FIPD,LIPD,FIAD,RV,Coup,DCC,EOM,DateOrigin,InputCheck=0,FindEOM,AnnivDatesOutput=BondAnalysis)
          if (BondAnalysis$Warnings$ZeroFlag==1) {
            Coup<-as.numeric(0)
            CpY<-as.numeric(1)
          }
          AccrInt<-DP.Output[[2]]$Accrued_Interest
          CF.values<-na.omit(BondAnalysis$PaySched$CoupPayments)
          CF.dates<-na.omit(BondAnalysis$PaySched$CoupDates)
          CF_List<-list(CF.values,CF.dates)
          Use.ClosedForm<-as.numeric(0)
          if (length(CF.values)>2) {
            if (length(which(!duplicated(CF.values[-c(1,length(CF.values))])==TRUE))==1) {
              Use.ClosedForm<-as.numeric(1)
            }
          }
          if (Calc.Method==0) {
            DCC_Orig<-DCC
            DCC<-2
            BondAnalysis<-suppressWarnings(AnnivDates(Em=Em,Mat=Mat,CpY=CpY,FIPD=FIPD,LIPD=LIPD,FIAD=FIAD,RV=RV,Coup=Coup,DCC=DCC,EOM=EOM,DateOrigin=DateOrigin,InputCheck=0,FindEOM=FindEOM,RegCF.equal=RegCF.equal))
          }
          # loading AnnivDates and expanding by one value to each side
          AD.set<-na.omit(BondAnalysis$DateVectors$AnnivDates)
          Refer<-BondAnalysis$Traits$Refer
          AtomVector_Refer<-as.numeric(unlist(strsplit(as.character(Refer),split = "-")))
          Atom1Refer<-AtomVector_Refer[1]
          Atom2Refer<-AtomVector_Refer[2]
          Atom3Refer<-AtomVector_Refer[3]
          # creating the anniversary date preceding AD1
          AtomVector_AD1<-as.numeric(unlist(strsplit(as.character(AD.set[1]),split = "-")))
          Atom1AD1<-AtomVector_AD1[1]
          Atom2AD1<-AtomVector_AD1[2]
          Atom3AD1<-AtomVector_AD1[3]
          PrevDate<-as.numeric(CppPrevDate(c(Atom1AD1,Atom2AD1,Atom3AD1,Atom1AD1,Atom2AD1,Atom3AD1,Atom1Refer,Atom2Refer,Atom3Refer,CpY,EOM)))
          PrevDate<-as.Date(paste(PrevDate[1],PrevDate[2],PrevDate[3],sep="-"))
          AD.set<-c(PrevDate,AD.set)
          AD.set<-sort(na.omit(AD.set[!duplicated(AD.set)]))
          # creating the anniversary date succeeding ADfin
          AtomVector_ADfin<-as.numeric(unlist(strsplit(as.character(AD.set[length(AD.set)]),split = "-")))
          Atom1ADfin<-AtomVector_ADfin[1]
          Atom2ADfin<-AtomVector_ADfin[2]
          Atom3ADfin<-AtomVector_ADfin[3]
          SuccDate<-as.numeric(CppSuccDate(c(Atom1ADfin,Atom2ADfin,Atom3ADfin,Atom1ADfin,Atom2ADfin,Atom3ADfin,Atom1Refer,Atom2Refer,Atom3Refer,CpY,EOM)))
          SuccDate<-as.Date(paste(SuccDate[1],SuccDate[2],SuccDate[3],sep="-"))
          AD.set<-c(AD.set,SuccDate)
          AD.set<-sort(na.omit(AD.set[!duplicated(AD.set)]))
          AD.indexes<-na.omit(BondAnalysis$DateVectors$AD_indexes)
          AD.indexes<-c((AD.indexes[1]-1),AD.indexes,(AD.indexes[length(AD.indexes)]+1))
          SD.set<-na.omit(BondAnalysis$DateVectors$RealDates)
          SD.indexes<-na.omit(BondAnalysis$DateVectors$RD_indexes)
          AD_List<-list(AD.set,AD.indexes)
          SD_List<-list(SD.set,SD.indexes)
          ### calculating tau
          AtomVector_Mat<-as.numeric(unlist(strsplit(as.character(Mat),split = "-")))
          Atom1Mat<-AtomVector_Mat[1]
          Atom2Mat<-AtomVector_Mat[2]
          Atom3Mat<-AtomVector_Mat[3]
          AtomVector_PCD.SETT<-as.numeric(unlist(strsplit(as.character(PCD(SETT,AD.set)),split = "-")))
          Atom1_PCD.SETT<-AtomVector_PCD.SETT[1]
          Atom2_PCD.SETT<-AtomVector_PCD.SETT[2]
          Atom3_PCD.SETT<-AtomVector_PCD.SETT[3]
          AtomVector_SETT<-as.numeric(unlist(strsplit(as.character(SETT),split = "-")))
          Atom1_SETT<-AtomVector_SETT[1]
          Atom2_SETT<-AtomVector_SETT[2]
          Atom3_SETT<-AtomVector_SETT[3]
          AtomVector_NCD.SETT<-as.numeric(unlist(strsplit(as.character(NCD(SETT,AD.set)),split = "-")))
          Atom1_NCD.SETT<-AtomVector_NCD.SETT[1]
          Atom2_NCD.SETT<-AtomVector_NCD.SETT[2]
          Atom3_NCD.SETT<-AtomVector_NCD.SETT[3]
          AtomVector_NCD.NCD.SETT<-as.numeric(unlist(strsplit(as.character(NCD(NCD(SETT,AD.set),AD.set)),split = "-")))
          Atom1_NCD.NCD.SETT<-AtomVector_NCD.NCD.SETT[1]
          Atom2_NCD.NCD.SETT<-AtomVector_NCD.NCD.SETT[2]
          Atom3_NCD.NCD.SETT<-AtomVector_NCD.NCD.SETT[3]
          if (is.element(DCC,c(1,3,5,6,8,10,11,12,15))) {
            tau_Num<-DIST(c(DCC,Atom1_PCD.SETT,Atom2_PCD.SETT,Atom3_PCD.SETT,Atom1_SETT,Atom2_SETT,Atom3_SETT))[2]
            tau_Den<-DIST(c(DCC,Atom1_PCD.SETT,Atom2_PCD.SETT,Atom3_PCD.SETT,Atom1_NCD.SETT,Atom2_NCD.SETT,Atom3_NCD.SETT))[2]
          }
          if (DCC==16) {
            NonBus.PCD.SETT<-length(which((NonBusDays.Brazil$Date>=PCD(SETT,AD.set))&(NonBusDays.Brazil$Date<SETT)))
            NonBus.PCD.NCD<-length(which((NonBusDays.Brazil$Date>=PCD(SETT,AD.set))&(NonBusDays.Brazil$Date<NCD(SETT,AD.set))))
            tau_Num<-DIST(c(DCC,Atom1_PCD.SETT,Atom2_PCD.SETT,Atom3_PCD.SETT,Atom1_SETT,Atom2_SETT,Atom3_SETT,NonBus.PCD.SETT))[2]
            tau_Den<-DIST(c(DCC,Atom1_PCD.SETT,Atom2_PCD.SETT,Atom3_PCD.SETT,Atom1_NCD.SETT,Atom2_NCD.SETT,Atom3_NCD.SETT,NonBus.PCD.NCD))[2]
          }
          if (DCC==2|DCC==14) {
            OrigDCC<-DCC
            DCC<-2
            tau_Num<-DIST(c(DCC,Atom1_PCD.SETT,Atom2_PCD.SETT,Atom3_PCD.SETT,Atom1_PCD.SETT,Atom2_PCD.SETT,Atom3_PCD.SETT,Atom1_NCD.SETT,Atom2_NCD.SETT,Atom3_NCD.SETT,
                            Atom1_PCD.SETT,Atom2_PCD.SETT,Atom3_PCD.SETT,Atom1_SETT,Atom2_SETT,Atom3_SETT,Atom1_NCD.SETT,Atom2_NCD.SETT,Atom3_NCD.SETT,
                            AD_List[[2]][which(AD_List[[1]]==PCD(SETT,AD.set))],
                            AD_List[[2]][which(AD_List[[1]]==NCD(PCD(SETT,AD.set),AD.set))],
                            CpY))[2]
            tau_Den<-DIST(c(DCC,Atom1_PCD.SETT,Atom2_PCD.SETT,Atom3_PCD.SETT,Atom1_PCD.SETT,Atom2_PCD.SETT,Atom3_PCD.SETT,Atom1_NCD.SETT,Atom2_NCD.SETT,Atom3_NCD.SETT,
                            Atom1_NCD.SETT,Atom2_NCD.SETT,Atom3_NCD.SETT,Atom1_NCD.SETT,Atom2_NCD.SETT,Atom3_NCD.SETT,Atom1_NCD.NCD.SETT,Atom2_NCD.NCD.SETT,Atom3_NCD.NCD.SETT,
                            AD_List[[2]][which(AD_List[[1]]==PCD(NCD(SETT,AD.set),AD.set))],
                            AD_List[[2]][which(AD_List[[1]]==NCD(PCD(SETT,AD.set),AD.set))],
                            CpY))[2]
            DCC<-OrigDCC
          }
          if (DCC==4) {
            tau_Num<-DIST(c(DCC,Atom1_PCD.SETT,Atom2_PCD.SETT,Atom3_PCD.SETT,Atom1_SETT,Atom2_SETT,Atom3_SETT,Atom1_NCD.SETT,CpY))[2]
            tau_Den<-DIST(c(DCC,Atom1_PCD.SETT,Atom2_PCD.SETT,Atom3_PCD.SETT,Atom1_NCD.SETT,Atom2_NCD.SETT,Atom3_NCD.SETT,Atom1_NCD.NCD.SETT,CpY))[2]
          }
          if (DCC==7) {
            tau_Num<-DIST(c(DCC,Atom1_PCD.SETT,Atom2_PCD.SETT,Atom3_PCD.SETT,Atom1_SETT,Atom2_SETT,Atom3_SETT,Atom1Mat,Atom2Mat,Atom3Mat))[2]
            tau_Den<-DIST(c(DCC,Atom1_PCD.SETT,Atom2_PCD.SETT,Atom3_PCD.SETT,Atom1_NCD.SETT,Atom2_NCD.SETT,Atom3_NCD.SETT,Atom1Mat,Atom2Mat,Atom3Mat))[2]
          }
          if (is.element(DCC,c(9,13))) {
            tau_Num<-DIST(c(DCC,Atom1_PCD.SETT,Atom2_PCD.SETT,Atom3_PCD.SETT,Atom1_SETT,Atom2_SETT,Atom3_SETT,EOM))[2]
            tau_Den<-DIST(c(DCC,Atom1_PCD.SETT,Atom2_PCD.SETT,Atom3_PCD.SETT,Atom1_NCD.SETT,Atom2_NCD.SETT,Atom3_NCD.SETT,EOM))[2]
          }
          tau<-tau_Num/tau_Den+AD_List[[2]][which(AD_List[[1]]==PCD(SETT,AD.set))]
          n<-SD.indexes[length(SD.indexes)-1]
          l<-BondAnalysis$Traits$LCPLength
          if ((!(is.na(LIPD)))&(SETT<LIPD)) {
            k<-SD_List[[2]][which(SD_List[[1]]==NCD(SETT,SD.set))]
            w<-k-tau
            eta<-n-k
            z<-l
            CN_tau<-CF_List[[1]][which(CF_List[[2]]==NCD(SETT,SD.set))]
          } else {
            k<-n
            w<-SD_List[[2]][which(SD_List[[1]]==NCD(SETT,SD.set))]-tau
            eta<-0
            z<-0
            CN_tau<-0
          }
          ###################>>>>>>                                               <<<<<<<##########################################################
          ##############>>>>>>        calculation of CP, ModDUR, MacDUR and Conv       <<<<<<<#####################################################
          ###################>>>>>>                                               <<<<<<<##########################################################
          CF_final<-CF.values[length(CF.values)]+RV
          if ((z==0)&(SimpleLastPeriod==TRUE)) { # simple interest
            a<-1+YtM*w/(100*CpY)
            DP<-CF_final/a
            CP<-DP-AccrInt
            ModDUR.inYears<-(-1)*(1/DP)*(-1)*CF_final*(a^(-2))*(w/CpY)
            MacDUR.inYears<-ModDUR.inYears*a
            Conv.inYears<-0.5*(1/DP)*2*CF_final*(a^(-3))*(w/CpY)^2
            ModDUR.inPeriods<-ModDUR.inYears*CpY
            MacDUR.inPeriods<-MacDUR.inYears*CpY
            Conv.inPeriods<-Conv.inYears*(CpY^2)
          } else {
            a<-1+YtM/(100*CpY)
            # if (round(a,8)==1) {
            #   CF.remain<-CF_List[[1]][which(CF_List[[2]]>SETT)]
            #   CF.remain<-c(CF.remain[-length(CF.remain)],CF.remain[length(CF.remain)]+RV)
            #   DP<-sum(CF.remain)
            #   CP<-DP-AccrInt
            # } else {
            #
            # }
            if ((Use.ClosedForm==1)&(!(round(a,8)==1))) {
              CF<-CF.values[2]
              DP<-(a^(-w))*(CN_tau+CF*(((a^eta)-1)/((a^eta)*(a-1)))+CF_final/(a^(eta+z)))
              CP<-DP-AccrInt
              ModDUR.inYears<-ModDUR(a,c(1,CN_tau,CF,CF_final,w,eta,z,CpY,DP))
              # ModDUR.inYears<-(-1)*(1/DP)*dm_MyPriceEqn(a,c(1,CN_tau,CF,CF_final,w,eta,z,CpY))
              MacDUR.inYears<-ModDUR.inYears*a
              Conv.inYears<-CONV(a,c(2,CN_tau,CF,CF_final,w,eta,z,CpY,DP))
              # Conv.inYears<-0.5*(1/DP)*dm_MyPriceEqn(a,c(2,CN_tau,CF,CF_final,w,eta,z,CpY))
              ModDUR.inPeriods<-ModDUR.inYears*CpY
              MacDUR.inPeriods<-MacDUR.inYears*CpY
              Conv.inPeriods<-Conv.inYears*(CpY^2)
            } else {
              if (z>0) {
                DiscPowerVector<-c((seq(0,eta,by=1)+w),(w+eta+z))
              } else {
                DiscPowerVector<-w
              }
              CF.remain<-CF_List[[1]][which(CF_List[[2]]>SETT)]
              CF.remain<-c(CF.remain[-length(CF.remain)],CF.remain[length(CF.remain)]+RV)
              # creating the price equation
              summands<-gsub(" ","",paste(CF.remain,"/(a^",DiscPowerVector,")"))
              PriceEqn<-paste(summands,collapse="+")
              PriceFunction.Standard<-function(a) {
                out<-eval(parse(text=PriceEqn))
                return(out)
              }
              # creating the first derivative of the price function
              d_PowerVector<-DiscPowerVector+1
              d_summands<-gsub(" ","",paste("(-1)*(",CF.remain,"*",DiscPowerVector,"/(a^",d_PowerVector,"))"))
              d_PriceEqn<-paste(d_summands,collapse="+")
              d_PriceFunction.Standard<-function(a) {
                out<-eval(parse(text=d_PriceEqn))
                return(out)
              }
              DP<-PriceFunction.Standard(a)
              CP<-DP-AccrInt
              # calculating the Modified Duration and the MacAulay Duration
              ModDUR.inPeriods<-(-1)*(1/DP)*d_PriceFunction.Standard(a)
              MacDUR.inPeriods<-ModDUR.inPeriods*a
              ModDUR.inYears<-ModDUR.inPeriods/CpY
              MacDUR.inYears<-MacDUR.inPeriods/CpY
              # calculating the Convexity
              d2_PowerVector<-d_PowerVector+1
              d2_summands<-gsub(" ","",paste("(",CF.remain,"*",DiscPowerVector,"*",d_PowerVector,")/(a^",d2_PowerVector,")"))
              Added.d2_summands<-paste(d2_summands,collapse="+")
              d2_PriceFunction.Standard <- function(a) {
                out<-eval(parse(text=Added.d2_summands))
                return(out)
              }
              Conv.inPeriods<-(1/DP)*0.5*d2_PriceFunction.Standard(a)
              Conv.inYears<-Conv.inPeriods/(CpY^2)
            }
          }
        }
      }
    }
  }
  # out<-list(tau=tau,n=n,k=k,w=w,eta=eta,z=z,CN_tau=CN_tau,ytm.p.a.=YtM,DP=DP,CF_final=CF_final,
  #           ModDUR.inPeriods=ModDUR.inPeriods,MacDUR.inPeriods=MacDUR.inPeriods,Conv.inPeriods=Conv.inPeriods,
  #           ModDUR.inYears=ModDUR.inYears,MacDUR.inYears=MacDUR.inYears,Conv.inYears=Conv.inYears,FIPD=FIPD,
  #           LIPD=LIPD)
  # out<-Added.d2_summands
  out<-list(CP=CP,AccrInt=AccrInt,DP=DP,ytm.p.a.=YtM,ModDUR.inYears=ModDUR.inYears,MacDUR.inYears=MacDUR.inYears,
            Conv.inYears=Conv.inYears,ModDUR.inPeriods=ModDUR.inPeriods,MacDUR.inPeriods=MacDUR.inPeriods,
            Conv.inPeriods=Conv.inPeriods,tau=tau)
  return(out)
}
