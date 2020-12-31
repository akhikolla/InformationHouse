
InputFormatCheck<-function(CP,YtM,SETT,Em,Mat,CpY,FIPD,LIPD,FIAD,RV,Coup,DCC,EOM,DateOrigin,StartDate,EndDate,YearNCP) {
  if (!missing(Em)) {
    if (!("Date"%in%class(Em))) {
      if ("numeric"%in%class(Em)) {
        Em<-as.Date(Em,origin=DateOrigin)
        warning(paste("The issue date (Em) is supplied as \"numeric\". It is converted to class \"Date\" using DateOrigin ",DateOrigin," and processed as Em =",Em,"."))
      } else {
        if (!(is.na(Em))) {
          if ("character"%in%class(Em)) {
            Em<-gsub("[ \\s]","",Em)
            HelpVector<-unlist(strsplit(Em,""))
            if (all(is.element(HelpVector,c("0":"9")))) {
              Em<-as.numeric(Em)
              Em<-as.Date(Em,origin=DateOrigin)
              warning(paste("The issue date (Em) is supplied as a number of class \"character\". It is converted to class \"Date\" using DateOrigin ",DateOrigin," and processed as Em =",Em,"."))
            } else {
              if (all(is.element(HelpVector[c(1:4,6,7,9,10)],c("0":"9")))) {
                if (HelpVector[5]=="-") {
                  if (HelpVector[8]=="-") {
                    Em<-as.Date(Em,"%Y-%m-%d")
                    if (is.na(Em)==FALSE) {
                      warning(paste("The issue date (Em) is supplied as a string of class \"character\" in the format \"yyyy-mm-dd\". It is converted to class \"Date\" using the command \"as.Date(Em,\"%Y-%m-%d\")\" and processed as Em =",Em,"."))
                    } else {
                      stop("The issue date (Em) cannot be processed!
                           Please make sure that it fits one of the following:
                           1. \"Date\" with format \"%Y-%m-%d\" or
                           2. \"numeric\" with the appropriate DateOrigin or
                           3. number of class \"character\" with the appropriate DateOrigin or
                           4. string of class \"character\" in the format \"yyyy-mm-dd\"
                           Note: If this date argument has class \"numeric\" or is a number of class \"character\",
                           origin can be set with the option DateOrigin (default is \"1970-01-01\").")
                    }
                  } else {
                    stop("The issue date (Em) cannot be processed!
                         Please make sure that it fits one of the following:
                         1. \"Date\" with format \"%Y-%m-%d\" or
                         2. \"numeric\" with the appropriate DateOrigin or
                         3. number of class \"character\" with the appropriate DateOrigin or
                         4. string of class \"character\" in the format \"yyyy-mm-dd\"
                         Note: If this date argument has class \"numeric\" or is a number of class \"character\",
                         origin can be set with the option DateOrigin (default is \"1970-01-01\").")
                  }
                } else {
                  stop("The issue date (Em) cannot be processed!
                       Please make sure that it fits one of the following:
                       1. \"Date\" with format \"%Y-%m-%d\" or
                       2. \"numeric\" with the appropriate DateOrigin or
                       3. number of class \"character\" with the appropriate DateOrigin or
                       4. string of class \"character\" in the format \"yyyy-mm-dd\"
                       Note: If this date argument has class \"numeric\" or is a number of class \"character\",
                       origin can be set with the option DateOrigin (default is \"1970-01-01\").")
                }
              } else {
                stop("The issue date (Em) cannot be processed!
                     Please make sure that it fits one of the following:
                     1. \"Date\" with format \"%Y-%m-%d\" or
                     2. \"numeric\" with the appropriate DateOrigin or
                     3. number of class \"character\" with the appropriate DateOrigin or
                     4. string of class \"character\" in the format \"yyyy-mm-dd\"
                     Note: If this date argument has class \"numeric\" or is a number of class \"character\",
                     origin can be set with the option DateOrigin (default is \"1970-01-01\").")
              }
            }
            rm(HelpVector)
          }
        }
      }
    }
  } else {
    Em=as.Date(NA)
  }
  if (!missing(Mat)) {
    if (!("Date"%in%class(Mat))) {
      if ("numeric"%in%class(Mat)) {
        Mat<-as.Date(Mat,origin=DateOrigin)
        warning(paste("The maturity date (Mat) is supplied as \"numeric\". It is converted to class \"Date\" using DateOrigin ",DateOrigin," and processed as Mat =",Mat,"."))
      } else {
        if (!(is.na(Mat))) {
          if ("character"%in%class(Mat)) {
            Mat<-gsub("[ \\s]","",Mat)
            HelpVector<-unlist(strsplit(Mat,""))
            if (all(is.element(HelpVector,c("0":"9")))) {
              Mat<-as.numeric(Mat)
              Mat<-as.Date(Mat,origin=DateOrigin)
              warning(paste("The maturity date (Mat) is supplied as a number of class \"character\". It is converted to class \"Date\" using DateOrigin ",DateOrigin," and processed as Mat =",Mat,"."))
            } else {
              if (all(is.element(HelpVector[c(1:4,6,7,9,10)],c("0":"9")))) {
                if (HelpVector[5]=="-") {
                  if (HelpVector[8]=="-") {
                    Mat<-as.Date(Mat,"%Y-%m-%d")
                    if (is.na(Mat)==FALSE) {
                      warning(paste("The maturity date (Mat) is supplied as a string of class \"character\" in the format \"yyyy-mm-dd\". It is converted to class \"Date\" using the command \"as.Date(Mat,\"%Y-%m-%d\")\" and processed as Mat =",Mat,"."))
                    } else {
                      stop("The maturity date (Mat) cannot be processed!
                           Please make sure that it fits one of the following:
                           1. \"Date\" with format \"%Y-%m-%d\" or
                           2. \"numeric\" with the appropriate DateOrigin or
                           3. number of class \"character\" with the appropriate DateOrigin or
                           4. string of class \"character\" in the format \"yyyy-mm-dd\"
                           Note: If this date argument has class \"numeric\" or is a number of class \"character\",
                           origin can be set with the option DateOrigin (default is \"1970-01-01\").")
                    }
                  } else {
                    stop("The maturity date (Mat) cannot be processed!
                         Please make sure that it fits one of the following:
                         1. \"Date\" with format \"%Y-%m-%d\" or
                         2. \"numeric\" with the appropriate DateOrigin or
                         3. number of class \"character\" with the appropriate DateOrigin or
                         4. string of class \"character\" in the format \"yyyy-mm-dd\"
                         Note: If this date argument has class \"numeric\" or is a number of class \"character\",
                         origin can be set with the option DateOrigin (default is \"1970-01-01\").")
                  }
                } else {
                  stop("The maturity date (Mat) cannot be processed!
                       Please make sure that it fits one of the following:
                       1. \"Date\" with format \"%Y-%m-%d\" or
                       2. \"numeric\" with the appropriate DateOrigin or
                       3. number of class \"character\" with the appropriate DateOrigin or
                       4. string of class \"character\" in the format \"yyyy-mm-dd\"
                       Note: If this date argument has class \"numeric\" or is a number of class \"character\",
                       origin can be set with the option DateOrigin (default is \"1970-01-01\").")
                }
              } else {
                stop("The maturity date (Mat) cannot be processed!
                     Please make sure that it fits one of the following:
                     1. \"Date\" with format \"%Y-%m-%d\" or
                     2. \"numeric\" with the appropriate DateOrigin or
                     3. number of class \"character\" with the appropriate DateOrigin or
                     4. string of class \"character\" in the format \"yyyy-mm-dd\"
                     Note: If this date argument has class \"numeric\" or is a number of class \"character\",
                     origin can be set with the option DateOrigin (default is \"1970-01-01\").")
              }
            }
            rm(HelpVector)
          }
        }
      }
    }
  } else {
    Mat=as.Date(NA)
  }
  if (!missing(FIPD)) {
    if (!("Date"%in%class(FIPD))) {
      if ("numeric"%in%class(FIPD)) {
        FIPD<-as.Date(FIPD,origin=DateOrigin)
        warning(paste("The first interest payment date (FIPD) is supplied as \"numeric\". It is converted to class \"Date\" using DateOrigin ",DateOrigin," and processed as FIPD =",FIPD,"."))
      } else {
        if (!(is.na(FIPD))) {
          if ("character"%in%class(FIPD)) {
            FIPD<-gsub("[ \\s]","",FIPD)
            HelpVector<-unlist(strsplit(FIPD,""))
            if (all(is.element(HelpVector,c("0":"9")))) {
              FIPD<-as.numeric(FIPD)
              FIPD<-as.Date(FIPD,origin=DateOrigin)
              warning(paste("The first interest payment date (FIPD) is supplied as a number of class \"character\". It is converted to class \"Date\" using DateOrigin ",DateOrigin," and processed as FIPD =",FIPD,"."))
            } else {
              if (all(is.element(HelpVector[c(1:4,6,7,9,10)],c("0":"9")))) {
                if (HelpVector[5]=="-") {
                  if (HelpVector[8]=="-") {
                    FIPD<-as.Date(FIPD,"%Y-%m-%d")
                    if (is.na(FIPD)==FALSE) {
                      warning(paste("The first interest payment date (FIPD) is supplied as a string of class \"character\" in the format \"yyyy-mm-dd\". It is converted to class \"Date\" using the command \"as.Date(FIPD,\"%Y-%m-%d\")\" and processed as FIPD =",FIPD,"."))
                    } else {
                      stop("The first interest payment date (FIPD) cannot be processed!
                           Please make sure that it fits one of the following:
                           1. \"Date\" with format \"%Y-%m-%d\" or
                           2. \"numeric\" with the appropriate DateOrigin or
                           3. number of class \"character\" with the appropriate DateOrigin or
                           4. string of class \"character\" in the format \"yyyy-mm-dd\"
                           Note: If this date argument has class \"numeric\" or is a number of class \"character\",
                           origin can be set with the option DateOrigin (default is \"1970-01-01\").")
                    }
                  } else {
                    stop("The first interest payment date (FIPD) cannot be processed!
                         Please make sure that it fits one of the following:
                         1. \"Date\" with format \"%Y-%m-%d\" or
                         2. \"numeric\" with the appropriate DateOrigin or
                         3. number of class \"character\" with the appropriate DateOrigin or
                         4. string of class \"character\" in the format \"yyyy-mm-dd\"
                         Note: If this date argument has class \"numeric\" or is a number of class \"character\",
                         origin can be set with the option DateOrigin (default is \"1970-01-01\").")
                  }
                } else {
                  stop("The first interest payment date (FIPD) cannot be processed!
                       Please make sure that it fits one of the following:
                       1. \"Date\" with format \"%Y-%m-%d\" or
                       2. \"numeric\" with the appropriate DateOrigin or
                       3. number of class \"character\" with the appropriate DateOrigin or
                       4. string of class \"character\" in the format \"yyyy-mm-dd\"
                       Note: If this date argument has class \"numeric\" or is a number of class \"character\",
                       origin can be set with the option DateOrigin (default is \"1970-01-01\").")
                }
              } else {
                stop("The first interest payment date (FIPD) cannot be processed!
                     Please make sure that it fits one of the following:
                     1. \"Date\" with format \"%Y-%m-%d\" or
                     2. \"numeric\" with the appropriate DateOrigin or
                     3. number of class \"character\" with the appropriate DateOrigin or
                     4. string of class \"character\" in the format \"yyyy-mm-dd\"
                     Note: If this date argument has class \"numeric\" or is a number of class \"character\",
                     origin can be set with the option DateOrigin (default is \"1970-01-01\").")
              }
            }
            rm(HelpVector)
          }
        }
      }
    }
  } else {
    FIPD=as.Date(NA)
  }
  if (!missing(LIPD)) {
    if (!("Date"%in%class(LIPD))) {
      if ("numeric"%in%class(LIPD)) {
        LIPD<-as.Date(LIPD,origin=DateOrigin)
        warning(paste("The last interest payment date (LIPD) is supplied as \"numeric\". It is converted to class \"Date\" using DateOrigin ",DateOrigin," and processed as LIPD =",LIPD,"."))
      } else {
        if (!(is.na(LIPD))) {
          if ("character"%in%class(LIPD)) {
            LIPD<-gsub("[ \\s]","",LIPD)
            HelpVector<-unlist(strsplit(LIPD,""))
            if (all(is.element(HelpVector,c("0":"9")))) {
              LIPD<-as.numeric(LIPD)
              LIPD<-as.Date(LIPD,origin=DateOrigin)
              warning(paste("The last interest payment date (LIPD) is supplied as a number of class \"character\". It is converted to class \"Date\" using DateOrigin ",DateOrigin," and processed as LIPD =",LIPD,"."))
            } else {
              if (all(is.element(HelpVector[c(1:4,6,7,9,10)],c("0":"9")))) {
                if (HelpVector[5]=="-") {
                  if (HelpVector[8]=="-") {
                    LIPD<-as.Date(LIPD,"%Y-%m-%d")
                    if (is.na(LIPD)==FALSE) {
                      warning(paste("The last interest payment date (LIPD) is supplied as a string of class \"character\" in the format \"yyyy-mm-dd\". It is converted to class \"Date\" using the command \"as.Date(LIPD,\"%Y-%m-%d\")\" and processed as LIPD =",LIPD,"."))
                    } else {
                      stop("The last interest payment date (LIPD) cannot be processed!
                           Please make sure that it fits one of the following:
                           1. \"Date\" with format \"%Y-%m-%d\" or
                           2. \"numeric\" with the appropriate DateOrigin or
                           3. number of class \"character\" with the appropriate DateOrigin or
                           4. string of class \"character\" in the format \"yyyy-mm-dd\"
                           Note: If this date argument has class \"numeric\" or is a number of class \"character\",
                           origin can be set with the option DateOrigin (default is \"1970-01-01\").")
                    }
                  } else {
                    stop("The last interest payment date (LIPD) cannot be processed!
                         Please make sure that it fits one of the following:
                         1. \"Date\" with format \"%Y-%m-%d\" or
                         2. \"numeric\" with the appropriate DateOrigin or
                         3. number of class \"character\" with the appropriate DateOrigin or
                         4. string of class \"character\" in the format \"yyyy-mm-dd\"
                         Note: If this date argument has class \"numeric\" or is a number of class \"character\",
                         origin can be set with the option DateOrigin (default is \"1970-01-01\").")
                  }
                } else {
                  stop("The last interest payment date (LIPD) cannot be processed!
                       Please make sure that it fits one of the following:
                       1. \"Date\" with format \"%Y-%m-%d\" or
                       2. \"numeric\" with the appropriate DateOrigin or
                       3. number of class \"character\" with the appropriate DateOrigin or
                       4. string of class \"character\" in the format \"yyyy-mm-dd\"
                       Note: If this date argument has class \"numeric\" or is a number of class \"character\",
                       origin can be set with the option DateOrigin (default is \"1970-01-01\").")
                }
              } else {
                stop("The last interest payment date (LIPD) cannot be processed!
                     Please make sure that it fits one of the following:
                     1. \"Date\" with format \"%Y-%m-%d\" or
                     2. \"numeric\" with the appropriate DateOrigin or
                     3. number of class \"character\" with the appropriate DateOrigin or
                     4. string of class \"character\" in the format \"yyyy-mm-dd\"
                     Note: If this date argument has class \"numeric\" or is a number of class \"character\",
                     origin can be set with the option DateOrigin (default is \"1970-01-01\").")
              }
            }
            rm(HelpVector)
          }
        }
      }
    }
  } else {
    LIPD=as.Date(NA)
  }
  if (!missing(FIAD)) {
    if (!("Date"%in%class(FIAD))) {
      if ("numeric"%in%class(FIAD)) {
        FIAD<-as.Date(FIAD,origin=DateOrigin)
        warning(paste("The first interest accrual date (FIAD) is supplied as \"numeric\". It is converted to class \"Date\" using DateOrigin ",DateOrigin," and processed as FIAD =",FIAD,"."))
      } else {
        if (!(is.na(FIAD))) {
          if ("character"%in%class(FIAD)) {
            FIAD<-gsub("[ \\s]","",FIAD)
            HelpVector<-unlist(strsplit(FIAD,""))
            if (all(is.element(HelpVector,c("0":"9")))) {
              FIAD<-as.numeric(FIAD)
              FIAD<-as.Date(FIAD,origin=DateOrigin)
              warning(paste("The first interest accrual date (FIAD) is supplied as a number of class \"character\". It is converted to class \"Date\" using DateOrigin ",DateOrigin," and processed as FIAD =",FIAD,"."))
            } else {
              if (all(is.element(HelpVector[c(1:4,6,7,9,10)],c("0":"9")))) {
                if (HelpVector[5]=="-") {
                  if (HelpVector[8]=="-") {
                    FIAD<-as.Date(FIAD,"%Y-%m-%d")
                    if (is.na(FIAD)==FALSE) {
                      warning(paste("The first interest accrual date (FIAD) is supplied as a string of class \"character\" in the format \"yyyy-mm-dd\". It is converted to class \"Date\" using the command \"as.Date(FIAD,\"%Y-%m-%d\")\" and processed as FIAD =",FIAD,"."))
                    } else {
                      stop("The first interest accrual date (FIAD) cannot be processed!
                           Please make sure that it fits one of the following:
                           1. \"Date\" with format \"%Y-%m-%d\" or
                           2. \"numeric\" with the appropriate DateOrigin or
                           3. number of class \"character\" with the appropriate DateOrigin or
                           4. string of class \"character\" in the format \"yyyy-mm-dd\"
                           Note: If this date argument has class \"numeric\" or is a number of class \"character\",
                           origin can be set with the option DateOrigin (default is \"1970-01-01\").")
                    }
                  } else {
                    stop("The first interest accrual date (FIAD) cannot be processed!
                         Please make sure that it fits one of the following:
                         1. \"Date\" with format \"%Y-%m-%d\" or
                         2. \"numeric\" with the appropriate DateOrigin or
                         3. number of class \"character\" with the appropriate DateOrigin or
                         4. string of class \"character\" in the format \"yyyy-mm-dd\"
                         Note: If this date argument has class \"numeric\" or is a number of class \"character\",
                         origin can be set with the option DateOrigin (default is \"1970-01-01\").")
                  }
                } else {
                  stop("The first interest accrual date (FIAD) cannot be processed!
                       Please make sure that it fits one of the following:
                       1. \"Date\" with format \"%Y-%m-%d\" or
                       2. \"numeric\" with the appropriate DateOrigin or
                       3. number of class \"character\" with the appropriate DateOrigin or
                       4. string of class \"character\" in the format \"yyyy-mm-dd\"
                       Note: If this date argument has class \"numeric\" or is a number of class \"character\",
                       origin can be set with the option DateOrigin (default is \"1970-01-01\").")
                }
              } else {
                stop("The first interest accrual date (FIAD) cannot be processed!
                     Please make sure that it fits one of the following:
                     1. \"Date\" with format \"%Y-%m-%d\" or
                     2. \"numeric\" with the appropriate DateOrigin or
                     3. number of class \"character\" with the appropriate DateOrigin or
                     4. string of class \"character\" in the format \"yyyy-mm-dd\"
                     Note: If this date argument has class \"numeric\" or is a number of class \"character\",
                     origin can be set with the option DateOrigin (default is \"1970-01-01\").")
              }
            }
            rm(HelpVector)
          }
        }
      }
    }
  } else {
    FIAD=as.Date(NA)
  }
  if (!missing(CpY)) {
    if (!("numeric"%in%class(CpY))) {
      if ("Date"%in%class(CpY)) {
        CpY<-as.numeric(CpY)
        warning(paste("The number of interest payments per year (CpY) is supplied as \"Date\".
                      Its conversion to class numeric results in CpY =",CpY,"."))
      } else {
        if (!(is.na(CpY))) {
          if ("character"%in%class(CpY)) {
            CpY<-gsub("[ \\s]","",CpY)
            HelpVector<-unlist(strsplit(CpY,""))
            if (all(is.element(HelpVector,c("0":"9")))) {
              CpY<-as.numeric(CpY)
              warning(paste("The number of interest payments per year (CpY) is supplied as a number of class \"character\".
                            Its conversion to class numeric results in CpY =",CpY,"."))
            } else {
              if (all(is.element(HelpVector[c(1:4,6,7,9,10)],c("0":"9")))) {
                if (HelpVector[5]=="-") {
                  if (HelpVector[8]=="-") {
                    CpY<-as.Date(CpY,"%Y-%m-%d")
                    CpY<-as.numeric(CpY)
                    if (!(is.na(CpY))) {
                      warning(paste("The number of interest payments per year (CpY) is supplied as a string of class \"character\"
                                    in the format \"yyyy-mm-dd\". It is converted to class \"Date\" using the command
                                    \"as.Date(CpY,\"%Y-%m-%d\")\" Its subsequent conversion to class numeric results in CpY =",CpY,"."))
                    } else {
                      stop("The number of interest payments per year (CpY) cannot be processed!
                           Please make sure that its class is either \"numeric\" or \"character\"!")
                    }
                  } else {
                    stop("The number of interest payments per year (CpY) cannot be processed!
                         Please make sure that its class is either \"numeric\" or \"character\"!")
                  }
                } else {
                  stop("The number of interest payments per year (CpY) cannot be processed!
                       Please make sure that its class is either \"numeric\" or \"character\"!")
                }
              } else {
                stop("The number of interest payments per year (CpY) cannot be processed!
                     Please make sure that its class is either \"numeric\" or \"character\"!")
              }
            }
            rm(HelpVector)
          }
        }
      }
    }
  } else {
    CpY=as.numeric(NA)
  }
  if (!missing(RV)) {
    if (!("numeric"%in%class(RV))) {
      if ("Date"%in%class(RV)) {
        RV<-as.numeric(RV)
        warning(paste("The redemption value (RV) is supplied as \"Date\".
                      Its conversion to class numeric results in RV =",RV,"."))
      } else {
        if (!(is.na(RV))) {
          if ("character"%in%class(RV)) {
            RV<-gsub("[ \\s]","",RV)
            HelpVector<-unlist(strsplit(RV,""))
            if ((all(is.element(HelpVector[-1],c(".","0":"9"))))&(is.element(HelpVector[1],c("-","0":"9")))) {
              RV<-as.numeric(RV)
              warning(paste("The redemption value (RV) is supplied as a number of class \"character\".
                            Its conversion to class numeric results in RV =",RV,"."))
            } else {
              if (all(is.element(HelpVector[c(1:4,6,7,9,10)],c(".","0":"9")))) {
                if (HelpVector[5]=="-") {
                  if (HelpVector[8]=="-") {
                    RV<-as.Date(RV,"%Y-%m-%d")
                    RV<-as.numeric(RV)
                    if (!(is.na(RV))) {
                      warning(paste("The redemption value (RV) is supplied as a string of class \"character\"
                                    in the format \"yyyy-mm-dd\". It is converted to class \"Date\" using the command
                                    \"as.Date(RV,\"%Y-%m-%d\")\" Its subsequent conversion to class numeric results in RV =",RV,"."))
                    } else {
                      stop("The redemption value (RV) cannot be processed!
                           Please make sure that its class is either \"numeric\" or \"character\"!")
                    }
                  } else {
                    stop("The redemption value (RV) cannot be processed!
                         Please make sure that its class is either \"numeric\" or \"character\"!")
                  }
                } else {
                  stop("The redemption value (RV) cannot be processed!
                       Please make sure that its class is either \"numeric\" or \"character\"!")
                }
              } else {
                stop("The redemption value (RV) cannot be processed!
                     Please make sure that its class is either \"numeric\" or \"character\"!")
              }
            }
            rm(HelpVector)
          }
        }
      }
    }
  } else {
    RV=as.numeric(NA)
  }
  if (!missing(Coup)) {
    if (!("numeric"%in%class(Coup))) {
      if ("Date"%in%class(Coup)) {
        Coup<-as.numeric(Coup)
        warning(paste("The nominal interest rate p.a. (Coup) is supplied as \"Date\".
                      Its conversion to class numeric results in Coup =",Coup,"."))
      } else {
        if (!(is.na(Coup))) {
          if ("character"%in%class(Coup)) {
            Coup<-gsub("[ \\s]","",Coup)
            HelpVector<-unlist(strsplit(Coup,""))
            if ((all(is.element(HelpVector[-1],c(".","0":"9"))))&(is.element(HelpVector[1],c("-","0":"9")))) {
              Coup<-as.numeric(Coup)
              warning(paste("The nominal interest rate p.a. (Coup) is supplied as a number of class \"character\".
                            Its conversion to class numeric results in Coup =",Coup,"."))
            } else {
              if (all(is.element(HelpVector[c(1:4,6,7,9,10)],c(".","0":"9")))) {
                if (HelpVector[5]=="-") {
                  if (HelpVector[8]=="-") {
                    Coup<-as.Date(Coup,"%Y-%m-%d")
                    Coup<-as.numeric(Coup)
                    if (!(is.na(Coup))) {
                      warning(paste("The nominal interest rate p.a. (Coup) is supplied as a string of class \"character\"
                                    in the format \"yyyy-mm-dd\". It is converted to class \"Date\" using the command
                                    \"as.Date(Coup,\"%Y-%m-%d\")\" Its subsequent conversion to class numeric results in Coup =",Coup,"."))
                    } else {
                      stop("The nominal interest rate p.a. (Coup) cannot be processed!
                           Please make sure that its class is either \"numeric\" or \"character\"!")
                    }
                  } else {
                    stop("The nominal interest rate p.a. (Coup) cannot be processed!
                         Please make sure that its class is either \"numeric\" or \"character\"!")
                  }
                } else {
                  stop("The nominal interest rate p.a. (Coup) cannot be processed!
                       Please make sure that its class is either \"numeric\" or \"character\"!")
                }
              } else {
                stop("The nominal interest rate p.a. (Coup) cannot be processed!
                     Please make sure that its class is either \"numeric\" or \"character\"!")
              }
            }
            rm(HelpVector)
          }
        }
      }
    }
  } else {
    Coup=as.numeric(NA)
  }
  if (!missing(DCC)) {
    if (!("numeric"%in%class(DCC))) {
      if ("Date"%in%class(DCC)) {
        DCC<-as.numeric(DCC)
        warning(paste("The day count convention identifier (DCC) is supplied as \"Date\".
                      Its conversion to class numeric results in DCC =",DCC,"."))
      } else {
        if (!(is.na(DCC))) {
          if ("character"%in%class(DCC)) {
            DCC<-gsub("[ \\s]","",DCC)
            HelpVector<-unlist(strsplit(DCC,""))
            if (all(is.element(HelpVector,c("0":"9")))) {
              DCC<-as.numeric(DCC)
              warning(paste("The day count convention identifier (DCC) is supplied as a number of class \"character\".
                            Its conversion to class numeric results in DCC =",DCC,"."))
            } else {
              if (all(is.element(HelpVector[c(1:4,6,7,9,10)],c("0":"9")))) {
                if (HelpVector[5]=="-") {
                  if (HelpVector[8]=="-") {
                    DCC<-as.Date(DCC,"%Y-%m-%d")
                    DCC<-as.numeric(DCC)
                    if (!(is.na(DCC))) {
                      warning(paste("The day count convention identifier (DCC) is supplied as a string of class \"character\"
                                    in the format \"yyyy-mm-dd\". It is converted to class \"Date\" using the command
                                    \"as.Date(DCC,\"%Y-%m-%d\")\" Its subsequent conversion to class numeric results in DCC =",DCC,"."))
                    } else {
                      stop("The day count convention identifier (DCC) cannot be processed!
                           Please make sure that its class is either \"numeric\" or \"character\"!")
                    }
                  } else {
                    stop("The day count convention identifier (DCC) cannot be processed!
                         Please make sure that its class is either \"numeric\" or \"character\"!")
                  }
                } else {
                  stop("The day count convention identifier (DCC) cannot be processed!
                       Please make sure that its class is either \"numeric\" or \"character\"!")
                }
              } else {
                stop("The day count convention identifier (DCC) cannot be processed!
                     Please make sure that its class is either \"numeric\" or \"character\"!")
              }
            }
            rm(HelpVector)
          }
        }
      }
    }
  } else {
    DCC=as.numeric(NA)
  }
  if (!(is.na(DCC))) {
    if (!(is.element(DCC,c(1:16)))) {
      stop("The day count convention identifier (DCC) cannot be processed!
           Please make sure that it is an element of the set {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16}!")
    }
  }
  if (!missing(EOM)) {
    if (!("numeric"%in%class(EOM))) {
      if ("Date"%in%class(EOM)) {
        EOM<-as.numeric(EOM)
        warning(paste("The End-of-Month rule boolean (EOM) is supplied as \"Date\".
                      Its conversion to class numeric results in EOM =",EOM,"."))
      } else {
        if (!(is.na(EOM))) {
          if ("character"%in%class(EOM)) {
            EOM<-gsub("[ \\s]","",EOM)
            HelpVector<-unlist(strsplit(EOM,""))
            if (all(is.element(HelpVector,c("0":"9")))) {
              EOM<-as.numeric(EOM)
              warning(paste("The End-of-Month rule boolean (EOM) is supplied as a number of class \"character\".
                            Its conversion to class numeric results in EOM =",EOM,"."))
            } else {
              if (all(is.element(HelpVector[c(1:4,6,7,9,10)],c("0":"9")))) {
                if (HelpVector[5]=="-") {
                  if (HelpVector[8]=="-") {
                    EOM<-as.Date(EOM,"%Y-%m-%d")
                    EOM<-as.numeric(EOM)
                    if (!(is.na(EOM))) {
                      warning(paste("The End-of-Month rule boolean (EOM) is supplied as a string of class \"character\"
                                    in the format \"yyyy-mm-dd\". It is converted to class \"Date\" using the command
                                    \"as.Date(EOM,\"%Y-%m-%d\")\" Its subsequent conversion to class numeric results in EOM =",EOM,"."))
                    } else {
                      stop("The End-of-Month rule boolean (EOM) cannot be processed!
                           Please make sure that its class is either \"numeric\" or \"character\"!")
                    }
                  } else {
                    stop("The End-of-Month rule boolean (EOM) cannot be processed!
                         Please make sure that its class is either \"numeric\" or \"character\"!")
                  }
                } else {
                  stop("The End-of-Month rule boolean (EOM) cannot be processed!
                       Please make sure that its class is either \"numeric\" or \"character\"!")
                }
              } else {
                stop("The End-of-Month rule boolean (EOM) cannot be processed!
                     Please make sure that its class is either \"numeric\" or \"character\"!")
              }
            }
            rm(HelpVector)
          }
        }
      }
    }
  } else {
    EOM=as.numeric(NA)
  }
  if (!(is.na(EOM))) {
    if (!(is.element(EOM,c(0,1)))) {
      stop("The End-of-Month rule boolean (EOM) cannot be processed!
           Please make sure that it is an element of the set {0,1}!")
    }
  }
  if (!missing(CP)) {
    if (!("numeric"%in%class(CP))) {
      if ("Date"%in%class(CP)) {
        CP<-as.numeric(CP)
        warning(paste("The clean price (CP) is supplied as \"Date\".
                      Its conversion to class numeric results in CP =",CP,"."))
      } else {
        if (!(is.na(CP))) {
          if ("character"%in%class(CP)) {
            CP<-gsub("[ \\s]","",CP)
            HelpVector<-unlist(strsplit(CP,""))
            if ((all(is.element(HelpVector[-1],c(".","0":"9"))))&(is.element(HelpVector[1],c("-","0":"9")))) {
              CP<-as.numeric(CP)
              warning(paste("The clean price (CP) is supplied as a number of class \"character\".
                            Its conversion to class numeric results in CP =",CP,"."))
            } else {
              if (all(is.element(HelpVector[c(1:4,6,7,9,10)],c(".","0":"9")))) {
                if (HelpVector[5]=="-") {
                  if (HelpVector[8]=="-") {
                    CP<-as.Date(CP,"%Y-%m-%d")
                    CP<-as.numeric(CP)
                    if (!(is.na(CP))) {
                      warning(paste("The clean price (CP) is supplied as a string of class \"character\"
                                    in the format \"yyyy-mm-dd\". It is converted to class \"Date\" using the command
                                    \"as.Date(CP,\"%Y-%m-%d\")\" Its subsequent conversion to class numeric results in CP =",CP,"."))
                    } else {
                      stop("The clean price (CP) cannot be processed!
                           Please make sure that its class is either \"numeric\" or \"character\"!")
                    }
                  } else {
                    stop("The clean price (CP) cannot be processed!
                         Please make sure that its class is either \"numeric\" or \"character\"!")
                  }
                } else {
                  stop("The clean price (CP) cannot be processed!
                       Please make sure that its class is either \"numeric\" or \"character\"!")
                }
              } else {
                stop("The clean price (CP) cannot be processed!
                     Please make sure that its class is either \"numeric\" or \"character\"!")
              }
            }
            rm(HelpVector)
          }
        }
      }
    }
  } else {
    CP=as.numeric(NA)
  }
  if (!missing(YtM)) {
    if (!("numeric"%in%class(YtM))) {
      if ("Date"%in%class(YtM)) {
        YtM<-as.numeric(YtM)
        warning(paste("The yield to maturity (YtM) is supplied as \"Date\".
                      Its conversion to class numeric results in YtM =",YtM,"."))
      } else {
        if (!(is.na(YtM))) {
          if ("character"%in%class(YtM)) {
            YtM<-gsub("[ \\s]","",YtM)
            HelpVector<-unlist(strsplit(YtM,""))
            if ((all(is.element(HelpVector[-1],c(".","0":"9"))))&(is.element(HelpVector[1],c("-","0":"9")))) {
              YtM<-as.numeric(YtM)
              warning(paste("The yield to maturity (YtM) is supplied as a number of class \"character\".
                            Its conversion to class numeric results in YtM =",YtM,"."))
            } else {
              if (all(is.element(HelpVector[c(1:4,6,7,9,10)],c(".","0":"9")))) {
                if (HelpVector[5]=="-") {
                  if (HelpVector[8]=="-") {
                    YtM<-as.Date(YtM,"%Y-%m-%d")
                    YtM<-as.numeric(YtM)
                    if (!(is.na(YtM))) {
                      warning(paste("The yield to maturity (YtM) is supplied as a string of class \"character\"
                                    in the format \"yyyy-mm-dd\". It is converted to class \"Date\" using the command
                                    \"as.Date(YtM,\"%Y-%m-%d\")\" Its subsequent conversion to class numeric results in YtM =",YtM,"."))
                    } else {
                      stop("The yield to maturity (YtM) cannot be processed!
                           Please make sure that its class is either \"numeric\" or \"character\"!")
                    }
                  } else {
                    stop("The yield to maturity (YtM) cannot be processed!
                         Please make sure that its class is either \"numeric\" or \"character\"!")
                  }
                } else {
                  stop("The yield to maturity (YtM) cannot be processed!
                       Please make sure that its class is either \"numeric\" or \"character\"!")
                }
              } else {
                stop("The yield to maturity (YtM) cannot be processed!
                     Please make sure that its class is either \"numeric\" or \"character\"!")
              }
            }
            rm(HelpVector)
          }
        }
      }
    }
  } else {
    YtM=as.numeric(NA)
  }
  if (!missing(SETT)) {
    if (!("Date"%in%class(SETT))) {
      if ("numeric"%in%class(SETT)) {
        SETT<-as.Date(SETT,origin=DateOrigin)
        warning(paste("The settlement date (SETT) is supplied as \"numeric\". It is converted to class \"Date\" using DateOrigin ",DateOrigin," and processed as SETT =",SETT,"."))
      } else {
        if (!(is.na(SETT))) {
          if ("character"%in%class(SETT)) {
            SETT<-gsub("[ \\s]","",SETT)
            HelpVector<-unlist(strsplit(SETT,""))
            if (all(is.element(HelpVector,c("0":"9")))) {
              SETT<-as.numeric(SETT)
              SETT<-as.Date(SETT,origin=DateOrigin)
              warning(paste("The settlement date (SETT) is supplied as a number of class \"character\". It is converted to class \"Date\" using DateOrigin ",DateOrigin," and processed as SETT =",SETT,"."))
            } else {
              if (all(is.element(HelpVector[c(1:4,6,7,9,10)],c("0":"9")))) {
                if (HelpVector[5]=="-") {
                  if (HelpVector[8]=="-") {
                    SETT<-as.Date(SETT,"%Y-%m-%d")
                    if (is.na(SETT)==FALSE) {
                      warning(paste("The settlement date (SETT) is supplied as a string of class \"character\" in the format \"yyyy-mm-dd\". It is converted to class \"Date\" using the command \"as.Date(SETT,\"%Y-%m-%d\")\" and processed as SETT =",SETT,"."))
                    } else {
                      stop("The settlement date (SETT) cannot be processed!
                           Please make sure that it fits one of the following:
                           1. \"Date\" with format \"%Y-%m-%d\" or
                           2. \"numeric\" with the appropriate DateOrigin or
                           3. number of class \"character\" with the appropriate DateOrigin or
                           4. string of class \"character\" in the format \"yyyy-mm-dd\"
                           Note: If this date argument has class \"numeric\" or is a number of class \"character\",
                           origin can be set with the option DateOrigin (default is \"1970-01-01\").")
                    }
                  } else {
                    stop("The settlement date (SETT) cannot be processed!
                         Please make sure that it fits one of the following:
                         1. \"Date\" with format \"%Y-%m-%d\" or
                         2. \"numeric\" with the appropriate DateOrigin or
                         3. number of class \"character\" with the appropriate DateOrigin or
                         4. string of class \"character\" in the format \"yyyy-mm-dd\"
                         Note: If this date argument has class \"numeric\" or is a number of class \"character\",
                         origin can be set with the option DateOrigin (default is \"1970-01-01\").")
                  }
                } else {
                  stop("The settlement date (SETT) cannot be processed!
                       Please make sure that it fits one of the following:
                       1. \"Date\" with format \"%Y-%m-%d\" or
                       2. \"numeric\" with the appropriate DateOrigin or
                       3. number of class \"character\" with the appropriate DateOrigin or
                       4. string of class \"character\" in the format \"yyyy-mm-dd\"
                       Note: If this date argument has class \"numeric\" or is a number of class \"character\",
                       origin can be set with the option DateOrigin (default is \"1970-01-01\").")
                }
              } else {
                stop("The settlement date (SETT) cannot be processed!
                     Please make sure that it fits one of the following:
                     1. \"Date\" with format \"%Y-%m-%d\" or
                     2. \"numeric\" with the appropriate DateOrigin or
                     3. number of class \"character\" with the appropriate DateOrigin or
                     4. string of class \"character\" in the format \"yyyy-mm-dd\"
                     Note: If this date argument has class \"numeric\" or is a number of class \"character\",
                     origin can be set with the option DateOrigin (default is \"1970-01-01\").")
              }
            }
            rm(HelpVector)
          }
        }
      }
    }
  } else {
    SETT=as.Date(NA)
  }
  if (!missing(StartDate)) {
    if (!("Date"%in%class(StartDate))) {
      if ("numeric"%in%class(StartDate)) {
        StartDate<-as.Date(StartDate,origin=DateOrigin)
        warning(paste("The starting date (StartDate) is supplied as \"numeric\". It is converted to class \"Date\" using DateOrigin ",DateOrigin," and processed as StartDate =",StartDate,"."))
      } else {
        if (!(is.na(StartDate))) {
          if ("character"%in%class(StartDate)) {
            StartDate<-gsub("[ \\s]","",StartDate)
            HelpVector<-unlist(strsplit(StartDate,""))
            if (all(is.element(HelpVector,c("0":"9")))) {
              StartDate<-as.numeric(StartDate)
              StartDate<-as.Date(StartDate,origin=DateOrigin)
              warning(paste("The starting date (StartDate) is supplied as a number of class \"character\". It is converted to class \"Date\" using DateOrigin ",DateOrigin," and processed as StartDate =",StartDate,"."))
            } else {
              if (all(is.element(HelpVector[c(1:4,6,7,9,10)],c("0":"9")))) {
                if (HelpVector[5]=="-") {
                  if (HelpVector[8]=="-") {
                    StartDate<-as.Date(StartDate,"%Y-%m-%d")
                    if (is.na(StartDate)==FALSE) {
                      warning(paste("The starting date (StartDate) is supplied as a string of class \"character\" in the format \"yyyy-mm-dd\". It is converted to class \"Date\" using the command \"as.Date(StartDate,\"%Y-%m-%d\")\" and processed as StartDate =",StartDate,"."))
                    } else {
                      stop("The starting date (StartDate) cannot be processed!
                           Please make sure that it fits one of the following:
                           1. \"Date\" with format \"%Y-%m-%d\" or
                           2. \"numeric\" with the appropriate DateOrigin or
                           3. number of class \"character\" with the appropriate DateOrigin or
                           4. string of class \"character\" in the format \"yyyy-mm-dd\"
                           Note: If this date argument has class \"numeric\" or is a number of class \"character\",
                           origin can be set with the option DateOrigin (default is \"1970-01-01\").")
                    }
                  } else {
                    stop("The starting date (StartDate) cannot be processed!
                         Please make sure that it fits one of the following:
                         1. \"Date\" with format \"%Y-%m-%d\" or
                         2. \"numeric\" with the appropriate DateOrigin or
                         3. number of class \"character\" with the appropriate DateOrigin or
                         4. string of class \"character\" in the format \"yyyy-mm-dd\"
                         Note: If this date argument has class \"numeric\" or is a number of class \"character\",
                         origin can be set with the option DateOrigin (default is \"1970-01-01\").")
                  }
                } else {
                  stop("The starting date (StartDate) cannot be processed!
                       Please make sure that it fits one of the following:
                       1. \"Date\" with format \"%Y-%m-%d\" or
                       2. \"numeric\" with the appropriate DateOrigin or
                       3. number of class \"character\" with the appropriate DateOrigin or
                       4. string of class \"character\" in the format \"yyyy-mm-dd\"
                       Note: If this date argument has class \"numeric\" or is a number of class \"character\",
                       origin can be set with the option DateOrigin (default is \"1970-01-01\").")
                }
              } else {
                stop("The starting date (StartDate) cannot be processed!
                     Please make sure that it fits one of the following:
                     1. \"Date\" with format \"%Y-%m-%d\" or
                     2. \"numeric\" with the appropriate DateOrigin or
                     3. number of class \"character\" with the appropriate DateOrigin or
                     4. string of class \"character\" in the format \"yyyy-mm-dd\"
                     Note: If this date argument has class \"numeric\" or is a number of class \"character\",
                     origin can be set with the option DateOrigin (default is \"1970-01-01\").")
              }
            }
            rm(HelpVector)
          }
        }
      }
    }
  } else {
    StartDate=as.Date(NA)
  }
  if (!missing(EndDate)) {
    if (!("Date"%in%class(EndDate))) {
      if ("numeric"%in%class(EndDate)) {
        EndDate<-as.Date(EndDate,origin=DateOrigin)
        warning(paste("The starting date (EndDate) is supplied as \"numeric\". It is converted to class \"Date\" using DateOrigin ",DateOrigin," and processed as EndDate =",EndDate,"."))
      } else {
        if (!(is.na(EndDate))) {
          if ("character"%in%class(EndDate)) {
            EndDate<-gsub("[ \\s]","",EndDate)
            HelpVector<-unlist(strsplit(EndDate,""))
            if (all(is.element(HelpVector,c("0":"9")))) {
              EndDate<-as.numeric(EndDate)
              EndDate<-as.Date(EndDate,origin=DateOrigin)
              warning(paste("The starting date (EndDate) is supplied as a number of class \"character\". It is converted to class \"Date\" using DateOrigin ",DateOrigin," and processed as EndDate =",EndDate,"."))
            } else {
              if (all(is.element(HelpVector[c(1:4,6,7,9,10)],c("0":"9")))) {
                if (HelpVector[5]=="-") {
                  if (HelpVector[8]=="-") {
                    EndDate<-as.Date(EndDate,"%Y-%m-%d")
                    if (is.na(EndDate)==FALSE) {
                      warning(paste("The starting date (EndDate) is supplied as a string of class \"character\" in the format \"yyyy-mm-dd\". It is converted to class \"Date\" using the command \"as.Date(EndDate,\"%Y-%m-%d\")\" and processed as EndDate =",EndDate,"."))
                    } else {
                      stop("The starting date (EndDate) cannot be processed!
                           Please make sure that it fits one of the following:
                           1. \"Date\" with format \"%Y-%m-%d\" or
                           2. \"numeric\" with the appropriate DateOrigin or
                           3. number of class \"character\" with the appropriate DateOrigin or
                           4. string of class \"character\" in the format \"yyyy-mm-dd\"
                           Note: If this date argument has class \"numeric\" or is a number of class \"character\",
                           origin can be set with the option DateOrigin (default is \"1970-01-01\").")
                    }
                  } else {
                    stop("The starting date (EndDate) cannot be processed!
                         Please make sure that it fits one of the following:
                         1. \"Date\" with format \"%Y-%m-%d\" or
                         2. \"numeric\" with the appropriate DateOrigin or
                         3. number of class \"character\" with the appropriate DateOrigin or
                         4. string of class \"character\" in the format \"yyyy-mm-dd\"
                         Note: If this date argument has class \"numeric\" or is a number of class \"character\",
                         origin can be set with the option DateOrigin (default is \"1970-01-01\").")
                  }
                } else {
                  stop("The starting date (EndDate) cannot be processed!
                       Please make sure that it fits one of the following:
                       1. \"Date\" with format \"%Y-%m-%d\" or
                       2. \"numeric\" with the appropriate DateOrigin or
                       3. number of class \"character\" with the appropriate DateOrigin or
                       4. string of class \"character\" in the format \"yyyy-mm-dd\"
                       Note: If this date argument has class \"numeric\" or is a number of class \"character\",
                       origin can be set with the option DateOrigin (default is \"1970-01-01\").")
                }
              } else {
                stop("The starting date (EndDate) cannot be processed!
                     Please make sure that it fits one of the following:
                     1. \"Date\" with format \"%Y-%m-%d\" or
                     2. \"numeric\" with the appropriate DateOrigin or
                     3. number of class \"character\" with the appropriate DateOrigin or
                     4. string of class \"character\" in the format \"yyyy-mm-dd\"
                     Note: If this date argument has class \"numeric\" or is a number of class \"character\",
                     origin can be set with the option DateOrigin (default is \"1970-01-01\").")
              }
            }
            rm(HelpVector)
          }
        }
      }
    }
  } else {
    EndDate=as.Date(NA)
  }
  if (!missing(YearNCP)) {
    if (!("numeric"%in%class(YearNCP))) {
      if ("Date"%in%class(YearNCP)) {
        YearNCP<-as.numeric(YearNCP)
        warning(paste("The year figure of next coupon date (YearNCP) is supplied as \"Date\".
                      Its conversion to class numeric results in YearNCP =",YearNCP,"."))
      } else {
        if (!(is.na(YearNCP))) {
          if ("character"%in%class(YearNCP)) {
            YearNCP<-gsub("[ \\s]","",YearNCP)
            HelpVector<-unlist(strsplit(YearNCP,""))
            if ((all(is.element(HelpVector[-1],c(".","0":"9"))))&(is.element(HelpVector[1],c("-","0":"9")))) {
              YearNCP<-as.numeric(YearNCP)
              warning(paste("The year figure of next coupon date (YearNCP) is supplied as a
                            number of class \"character\". Its conversion to class numeric
                            results in YearNCP =",YearNCP,"."))
            } else {
              if (all(is.element(HelpVector[c(1:4,6,7,9,10)],c(".","0":"9")))) {
                if (HelpVector[5]=="-") {
                  if (HelpVector[8]=="-") {
                    YearNCP<-as.Date(YearNCP,"%Y-%m-%d")
                    YearNCP<-as.numeric(YearNCP)
                    if (!(is.na(YearNCP))) {
                      warning(paste("The year figure of next coupon date (YearNCP) is supplied as a string
                                    of class \"character\" in the format \"yyyy-mm-dd\". It is converted to
                                    class \"Date\" using the command \"as.Date(YearNCP,\"%Y-%m-%d\")\" Its
                                    subsequent conversion to class numeric results in YearNCP =",YearNCP,"."))
                    } else {
                      stop("The year figure of next coupon date (YearNCP) cannot be processed!
                           Please make sure that its class is either \"numeric\" or \"character\"!")
                    }
                  } else {
                    stop("The year figure of next coupon date (YearNCP) cannot be processed!
                         Please make sure that its class is either \"numeric\" or \"character\"!")
                  }
                } else {
                  stop("The year figure of next coupon date (YearNCP) cannot be processed!
                       Please make sure that its class is either \"numeric\" or \"character\"!")
                }
              } else {
                stop("The year figure of next coupon date (YearNCP) cannot be processed!
                     Please make sure that its class is either \"numeric\" or \"character\"!")
              }
            }
            rm(HelpVector)
          }
        }
      }
    }
  } else {
    YearNCP=as.numeric(NA)
  }
  if (!(is.na(YearNCP))) {
    if (abs(YearNCP-round(YearNCP))>.Machine$double.eps^0.5) {
      stop("The year figure of next coupon date (YearNCP) cannot be processed!
           Please make sure that it is a year figure!")
    }
  }
  CheckedInput<-list(CP=CP,YtM=YtM,SETT=SETT,Em=Em,Mat=Mat,CpY=CpY,FIPD=FIPD,LIPD=LIPD,FIAD=FIAD,RV=RV,Coup=Coup,DCC=DCC,
                     EOM=EOM,DateOrigin=DateOrigin,StartDate=StartDate,EndDate=EndDate,YearNCP=YearNCP)
  return(CheckedInput)
}


PCD<-function(Date,DateVector) {
  if (!(("Date"%in%class(Date))|("numeric"%in%class(Date)))) {
    stop("The scalar cannot be processed! Please make sure that it fits one of the following: 1. \"Date\" with format \"%Y-%m-%d\" or 2. \"numeric\".")
  } else {
    if (!(("Date"%in%class(DateVector))|("numeric"%in%class(DateVector)))) {
      stop("The vector cannot be processed! Please make sure that it fits one of the following: 1. \"Date\" with format \"%Y-%m-%d\" or 2. \"numeric\".")
    } else {
      DateVector<-sort(na.omit(DateVector[!duplicated(DateVector)]))
      if (length(DateVector)<2) {
        stop("The vector must contain at least two different non-NA elements!")
      } else {
        if ((Date<DateVector[1])|(DateVector[length(DateVector)]<Date)) {
          PrevDate<-NA
          warning("The entered scalar is outside the vector boundaries. NA created!")
        } else {
          DateVector<-append(DateVector,Date)
          DateVector<-sort(DateVector)
          if (anyDuplicated(DateVector)!=0) {
            PrevDate<-DateVector[anyDuplicated(DateVector)]
          } else {
            PrevDate<-DateVector[which(DateVector<Date)[length(DateVector[which(DateVector<Date)])]]
          }
        }
        return(PrevDate)
      }
    }
  }
}


NCD<-function(Date,DateVector) {
  if (!(("Date"%in%class(Date))|("numeric"%in%class(Date)))) {
    stop("The scalar cannot be processed! Please make sure that it fits one of the following: 1. \"Date\" with format \"%Y-%m-%d\" or 2. \"numeric\".")
  } else {
    if (!(("Date"%in%class(DateVector))|("numeric"%in%class(DateVector)))) {
      stop("The vector cannot be processed! Please make sure that it fits one of the following: 1. \"Date\" with format \"%Y-%m-%d\" or 2. \"numeric\".")
    } else {
      DateVector<-sort(na.omit(DateVector[!duplicated(DateVector)]))
      if (length(DateVector)<2) {
        stop("The vector must contain at least two different non-NA elements!")
      } else {
        if ((Date<DateVector[1])|(!(DateVector[length(DateVector)]>Date))) {
          NextDate<-NA
          warning("The entered scalar is outside the vector boundaries. NA created!")
        } else {
          DateVector<-append(DateVector,Date)
          DateVector<-sort(DateVector)
          NextDate<-DateVector[which(DateVector>Date)[1]]
        }
        return(NextDate)
      }
    }
  }
}


NewtonRaphson<-function(f,df,startValue,Precision=.Machine$double.eps^0.25) {
  i<-0
  v<-startValue
  vNext<-v-f(v)/df(v)
  path<-c(v,vNext)
  path.break<-as.numeric(0)
  f.value<-f(vNext)
  while ((!(is.nan(f.value)))&(abs(f.value)>Precision)) {
    i<-i+1
    v<-vNext
    vNext<-v-f(v)/df(v)
    path<-append(path,vNext)
    f.value<-f(vNext)
    # if (round(abs(path[length(path)]-path[length(path)-1])-
    #           abs(path[length(path)-1]-path[length(path)-2]))<Precision) {
    #   path.break<-as.numeric(1)
    #   break
    # }
  }
  # if ((1+0.0000001>vNext)&(1-0.0000001<vNext)) {
  #   path.break<-as.numeric(1)
  # }
  out<-list(root=vNext,N.Iter=i,f.value=f.value,path.break=path.break)
  return(out)
}


dm_MyPriceEqn<-function(a,x) {
  m<-x[1]
  CNtau<-x[2]
  CF<-x[3]
  CF_final<-x[4]
  w<-x[5]
  eta<-x[6]
  z<-x[7]
  CpY<-x[8]
  Factor1<-((-1)^m)/((a^(w+m))*(CpY^(m)))
  Summand1<-CNtau*gamma(w+m)/(gamma(w))
  InnerFactor<-CF*factorial(m)/(a-1)
  mSequence<-seq(0,m)
  InnerSummands<-(1/(factorial(m-mSequence)))*((a/(a-1))^mSequence)*
    (((gamma(w+m-mSequence))/(gamma(w)))-
       (gamma(w+eta+m-mSequence))/((a^eta)*(gamma(w+eta))))
  InnerSum<-sum(InnerSummands)
  LastSummand<-(CF_final*gamma(w+eta+z+m))/((a^(eta+z))*(gamma(w+eta+z)))
  out<-Factor1*(Summand1+InnerFactor*InnerSum+LastSummand)
  return(out)
}


ModDUR<-function(a,x) {
  m<-x[1]
  CNtau<-x[2]
  CF<-x[3]
  CF_final<-x[4]
  w<-x[5]
  eta<-x[6]
  z<-x[7]
  CpY<-x[8]
  DP<-x[9]
  Factor1<-(1/(DP*CpY*a^(w+1)))
  Summand1<-CNtau*w
  InnerFactor<-CF*((a^eta)-1)/((a^eta)*(a-1))
  InnerSum<-w-(eta/((a^eta)-1))+(a/(a-1))
  LastSummand<-(CF_final*(w+eta+z))/(a^(eta+z))
  out<-Factor1*(Summand1+InnerFactor*InnerSum+LastSummand)
  return(out)
}


CONV<-function(a,x) {
  m<-x[1]
  CNtau<-x[2]
  CF<-x[3]
  CF_final<-x[4]
  w<-x[5]
  eta<-x[6]
  z<-x[7]
  CpY<-x[8]
  DP<-x[9]
  Factor1<-(1/(2*(a^w)*DP))*((1/(CpY*a))^2)
  Summand1<-CNtau*w*(w+1)
  InnerFactor<-CF*((a^eta)-1)/((a^eta)*(a-1))
  InnerSum<-w*(w+1)-
    ((eta*(2*w+eta+1))/((a^eta)-1))+
    ((2*a)/(a-1))*
    (w-eta/(a^eta-1)+(a/(a-1)))
  LastSummand<-(CF_final*(w+eta+z)*(w+eta+z+1))/(a^(eta+z))
  out<-Factor1*(Summand1+InnerFactor*InnerSum+LastSummand)
  return(out)
}

