#include <math.h>
#include <Rcpp.h>
#include <numeric>
#include <cmath>
using namespace Rcpp;
using namespace std;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
int leap (int x) {
  int year = x;
  int rest400 = year % 400;
  int rest100 = year % 100;
  int rest4 = year % 4;
  int out;
  if (rest400==0) {
    out = 1;
  } else {
    if (rest100==0) {
      out = 0;
    } else {
      if (rest4==0) {
        out = 1;
      } else {
        out = 0;
      }
    }
  }
  return out;
}



// [[Rcpp::export]]
int LDM (IntegerVector x) {
  int Y = x[0];
  int M = x[1];
  int D = x[2];
  int LDM;
  if (M==1) {
    if (D==31) {
      LDM = 1;
    } else {
      LDM = 0;
    }
  } else {
    if (M==2) {
      int LEAP = leap(Y);
      if (LEAP==1) {
        if (D==29) {
          LDM = 1;
        } else {
          LDM = 0;
        }
      } else {
        if (D==28) {
          LDM = 1;
        } else {
          LDM = 0;
        }
      }
    } else {
      if (M==3) {
        if (D==31) {
          LDM = 1;
        } else {
          LDM = 0;
        }
      } else {
        if (M==4) {
          if (D==30) {
            LDM = 1;
          } else {
            LDM = 0;
          }
        } else {
          if (M==5) {
            if (D==31) {
              LDM = 1;
            } else {
              LDM = 0;
            }
          } else {
            if (M==6) {
              if (D==30) {
                LDM = 1;
              } else {
                LDM = 0;
              }
            } else {
              if (M==7) {
                if (D==31) {
                  LDM = 1;
                } else {
                  LDM = 0;
                }
              } else {
                if (M==8) {
                  if (D==31) {
                    LDM = 1;
                  } else {
                    LDM = 0;
                  }
                } else {
                  if (M==9) {
                    if (D==30) {
                      LDM = 1;
                    } else {
                      LDM = 0;
                    }
                  } else {
                    if (M==10) {
                      if (D==31) {
                        LDM = 1;
                      } else {
                        LDM = 0;
                      }
                    } else {
                      if (M==11) {
                        if (D==30) {
                          LDM = 1;
                        } else {
                          LDM = 0;
                        }
                      } else { // i.e. if M==12
                        if (D==31) {
                          LDM = 1;
                        } else {
                          LDM = 0;
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
  return LDM;
}



// [[Rcpp::export]]
int DaysInMonth (IntegerVector x) {
  int Y = x[0];
  int M = x[1];
  int LEAP = leap(Y);
  int out;
  if (M==1) {
    out = 31;
  } else {
    if (M==2) {
      if (LEAP==1) {
        out = 29;
      } else {
        out = 28;
      }
    } else {
      if (M==3) {
        out = 31;
      } else {
        if (M==4) {
          out = 30;
        } else {
          if (M==5) {
            out = 31;
          } else {
            if (M==6) {
              out = 30;
            } else {
              if (M==7) {
                out = 31;
              } else {
                if (M==8) {
                  out = 31;
                } else {
                  if (M==9) {
                    out = 30;
                  } else {
                    if (M==10) {
                      out = 31;
                    } else {
                      if (M==11) {
                        out = 30;
                      } else { // i.e. if M==12
                        out = 31;
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
  return out;
}



// [[Rcpp::export]]
int DaysInYear(int x) {
  int LEAP = leap(x);
  int out;
  if (LEAP==1) {
    out = 366;
  } else {
    out = 365;
  }
  return out;
}



// [[Rcpp::export]]
int DayDiff(IntegerVector x) {
  int Y1 = x[0];
  int M1 = x[1];
  int D1 = x[2];
  int Y2 = x[3];
  int M2 = x[4];
  int D2 = x[5];
  int DateSwitch = 0;
  if (x[0]>x[3]) {
    Y2 = x[0];
    M2 = x[1];
    D2 = x[2];
    Y1 = x[3];
    M1 = x[4];
    D1 = x[5];
    DateSwitch = 1;
  } else {
    if (x[0]==x[3]) {
      if (x[1]>x[4]) {
        Y2 = x[0];
        M2 = x[1];
        D2 = x[2];
        Y1 = x[3];
        M1 = x[4];
        D1 = x[5];
        DateSwitch = 1;
      } else {
        if (x[1]==x[4]) {
          if (x[2]>x[5]) {
            Y2 = x[0];
            M2 = x[1];
            D2 = x[2];
            Y1 = x[3];
            M1 = x[4];
            D1 = x[5];
            DateSwitch = 1;
          }
        }
      }
    }
  }
  IntegerVector Date1(3);
  IntegerVector Date2(3);
  int values1[] = {Y1,M1,D1};
  int values2[] = {Y2,M2,D2};
  Date1.assign(values1,values1+3);
  Date2.assign(values2,values2+3);
  int DaysInM2 = DaysInMonth(Date2);
  int DaysInY1 = DaysInYear(Y1);
  int DaysInY2 = DaysInYear(Y2);
  int Ydiff = Y2 - Y1;
  int Mdiff;
  int out;
  if (Ydiff==0) {
    Mdiff = M2 - M1;
    if (Mdiff==0) {
      out = D2 - D1;
    } else {
      IntegerVector MVector(Mdiff+1);
      IntegerVector CalcVector(2);
      CalcVector[0] = Y1;
      for(int i = 0; i < Mdiff + 1; ++i) {
        CalcVector[1] = M1 + i;
        MVector[i] = DaysInMonth(CalcVector);
      }
      int Msum = accumulate(MVector.begin(), MVector.end(), 0);
      out = Msum - D1 + D2 - DaysInM2;
    }
  } else {
    IntegerVector YVector(Ydiff+1);
    for(int i = 0; i < Ydiff + 1; ++i) {
      YVector[i] = DaysInYear(Y1 + i);
    }
    int Ysum = accumulate(YVector.begin(), YVector.end(), 0);
    int StartYMdiff = 12 - M1;
    IntegerVector StartYMVector(StartYMdiff+1);
    IntegerVector CalcVector(2);
    CalcVector[0] = Y1;
    for(int i = 0; i < StartYMdiff + 1; ++i) {
      CalcVector[1] = M1 + i;
      StartYMVector[i] = DaysInMonth(CalcVector);
    }
    int StartYMsum = accumulate(StartYMVector.begin(), StartYMVector.end(), 0);
    IntegerVector EndYMVector(M2);
    IntegerVector EndCalcVector(2);
    EndCalcVector[0] = Y2;
    for(int i = 0; i < M2; ++i) {
      EndCalcVector[1] = i + 1;
      EndYMVector[i] = DaysInMonth(EndCalcVector);
    }
    int EndYMsum = accumulate(EndYMVector.begin(), EndYMVector.end(), 0);
    out = Ysum - DaysInY1 - DaysInY2 + StartYMsum - D1  + EndYMsum - DaysInM2 + D2;
  }
  if (DateSwitch==1) {
    out = (-1)*out;
  }
  return out;
}



// [[Rcpp::export]]
IntegerVector Date_LDM(IntegerVector x) {
  IntegerVector out = x;
  int out_day = DaysInMonth(x);
  out(2) = out_day;
  return out;
}



// [[Rcpp::export]]
double sumC(NumericVector x) {
  int n = x.size();
  double total = 0;
  for(int i = 0; i < n; ++i) {
    total += x[i];
  }
  return total;
}



// [[Rcpp::export]]
int FirstMatch(IntegerVector x) {
  // first element in x is the element to be matched
  // all other elements of x are the vector to search
  // function returns the position of the first occurance of element in vector
  // if element is not in vector, function returns (vector.size+1)
  int Element = x[0];
  IntegerVector Menge(x.begin()+1,x.end());
  int Position = find(Menge.begin(), Menge.end(), Element) - Menge.begin() + 1;
  return Position;
}



// [[Rcpp::export]]
int LeapDayInside(IntegerVector x) {
  int Y1 = x[0];
  int M1 = x[1];
  int D1 = x[2];
  int Y2 = x[3];
  int M2 = x[4];
  int D2 = x[5];
  int out;
  int Ydiff = Y2 - Y1;
  if (Ydiff==0) {
    if (leap(Y1)==1) {
      int Date1 = 10000*Y1 + 100*M1 + D1;
      int Date2 = 10000*Y2 + 100*M2 + D2;
      int LeapDay = 10000*Y1 + 229;
      if ((Date1<=LeapDay)&(LeapDay<Date2)) {
        out = 1;
      } else {
        out = 0;
      }
    } else {
      out = 0;
    }
  } else {
    int lengthLeapCheckVec = Ydiff + 1;
    IntegerVector LeapCheckVec(lengthLeapCheckVec);
    for(int i = 0; i < lengthLeapCheckVec; ++i) {
      LeapCheckVec[i] = leap(Y1 + i);
    }
    if (sumC((NumericVector)LeapCheckVec)>2) {
      out = 1;
    } else {
      if (sumC((NumericVector)LeapCheckVec)==0) {
        out = 0;
      } else {
        if (sumC((NumericVector)LeapCheckVec)==1) {
          int WhichIsLeap = find(LeapCheckVec.begin(), LeapCheckVec.end(), 1) - LeapCheckVec.begin() + 1;
          if (WhichIsLeap==1) {
            int Date1 = 10000*Y1 + 100*M1 + D1;
            int Date2 = 10000*Y2 + 100*M2 + D2;
            int LeapDay = 10000*Y1 + 229;
            if ((Date1<=LeapDay)&(LeapDay<Date2)) {
              out = 1;
            } else {
              out = 0;
            }
          } else {
            if (WhichIsLeap==lengthLeapCheckVec) {
              int Date1 = 10000*Y1 + 100*M1 + D1;
              int Date2 = 10000*Y2 + 100*M2 + D2;
              int LeapDay = 10000*Y2 + 229;
              if ((Date1<=LeapDay)&(LeapDay<Date2)) {
                out = 1;
              } else {
                out = 0;
              }
            } else {
              out = 1;
            }
          }
        } else {
          if (LeapCheckVec.size()==5) {
            int Date1 = 10000*Y1 + 100*M1 + D1;
            int Date2 = 10000*Y2 + 100*M2 + D2;
            int LeapDay1 = 10000*Y1 + 229;
            int LeapDay2 = 10000*Y2 + 229;
            if ((Date1<=LeapDay1)|(LeapDay2<Date2)) {
              out = 1;
            } else {
              out = 0;
            }
          } else {
            out = 0;
          }
        }
      }
    }
  }
  return out;
}



// [[Rcpp::export]]
NumericVector DIST(NumericVector x) {
  // for DCC = {1,3,5,6,8,10,11,12,15} x is a vector of 7 integers:
  // c(DCC,Y1,M1,D1,Y2,M2,D2) with
  // Y1-M1-D1 = t_a ; Y2-M2-D2 = t_b
  //
  // for DCC = {2,14} x is a vector of 22 integers:
  // c(DCC,Y1,M1,D1,Y2,M2,D2,Y3,M3,D3,Y4,M4,D4,Y5,M5,D5,Y6,M6,D6,P,N,CpY) with
  // Y1-M1-D1 = PCD(t_a,AD) ; Y2-M2-D2 = t_a ; Y3-M3-D3 = NCD(t_a,AD)
  // Y4-M4-D4 = PCD(t_b,AD) ; Y5-M5-D5 = t_b ; Y6-M6-D6 = NCD(t_b,AD)
  // P = P(t_b,AD) ; N = N(t_a,AD)
  //
  // for DCC = 4 x is a vector of 9 integers:
  // c(DCC,Y1,M1,D1,Y2,M2,D2,Y3,CpY) with
  // Y1-M1-D1 = t_a ; Y2-M2-D2 = t_b;
  // Y3 = Year figure of the next coupon payment date after t_b
  //
  // for DCC = 7 x is a vector of 10 integers:
  // c(DCC,Y1,M1,D1,Y2,M2,D2,Y3,M3,D3) with
  // Y1-M1-D1 = t_a ; Y2-M2-D2 = t_b ; Y3-M3-D3 = t_M
  //
  // for DCC = {9,13} x is a vector of 8 integers:
  // c(DCC,Y1,M1,D1,Y2,M2,D2,EOM) with
  // Y1-M1-D1 = t_a ; Y2-M2-D2 = t_b
  //
  int DCC = x[0];
  int Y1 = x[1];
  int M1 = x[2];
  int D1 = x[3];
  int Y2 = x[4];
  int M2 = x[5];
  int D2 = x[6];
  double NAccr;
  double DIST_value;
  NumericVector out;
  if (DCC==1) { // Act/Act ISDA; DS-Var-BAB = 6
    // in this case R should pass a vector of 7 integers:
    // c(DCC,Y1,M1,D1,Y2,M2,D2) with
    // Y1-M1-D1 = t_a ; Y2-M2-D2 = t_b
    IntegerVector Diff12(6);
    int valuesDiff12[] = {Y1,M1,D1,Y2,M2,D2};
    Diff12.assign(valuesDiff12,valuesDiff12+6);
    int Ydiff21 = Y2 - Y1;
    int DaysInY1 = DaysInYear(Y1);
    if (Ydiff21==0) {
      NAccr = DayDiff(Diff12);
      DIST_value = (double)NAccr / (double)DaysInY1;
    } else {
      int lengthYVect21 = Ydiff21 + 1;
      IntegerVector YVector(lengthYVect21);
      for(int i = 0; i < lengthYVect21; ++i) {
        YVector[i] = DaysInYear(Y1 + i);
      }
      int valuesDaysY1[] = {Y1,M1,D1,Y1+1,1,1};
      IntegerVector DaysY1Vector(6);
      DaysY1Vector.assign(valuesDaysY1,valuesDaysY1+6);
      int DaysY1 = DayDiff(DaysY1Vector);
      int valuesDaysY2[] = {Y2,1,1,Y2,M2,D2};
      IntegerVector DaysY2Vector(6);
      DaysY2Vector.assign(valuesDaysY2,valuesDaysY2+6);
      int DaysY2 = DayDiff(DaysY2Vector);
      IntegerVector AccrDayVec(lengthYVect21);
      AccrDayVec[0] = DaysY1;
      for(int i = 1; i < lengthYVect21-1; ++i) {
        AccrDayVec[i] = YVector[i];
      }
      AccrDayVec[lengthYVect21-1] = DaysY2;
      NAccr = accumulate(AccrDayVec.begin(), AccrDayVec.end(), 0);
      NumericVector AccrFacVec(lengthYVect21);
      for (int i = 0; i < lengthYVect21; ++i) {
        double Num = AccrDayVec[i];
        double Den = YVector[i];
        AccrFacVec[i] = Num / Den;
      }
      DIST_value = sumC(AccrFacVec);
    }
  } else {
    if (DCC==2) { // Act/Act (ICMA); Act/Act (ISMA); DS-Var-BAB = 5
      // in this case R should pass a vector of 22 integers:
      // c(DCC,Y1,M1,D1,Y2,M2,D2,Y3,M3,D3,Y4,M4,D4,Y5,M5,D5,Y6,M6,D6,P,N,CpY) with
      // Y1-M1-D1 = PCD(t_a,AD) ; Y2-M2-D2 = t_a ; Y3-M3-D3 = NCD(t_a,AD)
      // Y4-M4-D4 = PCD(t_b,AD) ; Y5-M5-D5 = t_b ; Y6-M6-D6 = NCD(t_b,AD)
      // P = P(t_b,AD) ; N = N(t_a,AD)
      int Y3 = x[7];
      int M3 = x[8];
      int D3 = x[9];
      int Y4 = x[10];
      int M4 = x[11];
      int D4 = x[12];
      int Y5 = x[13];
      int M5 = x[14];
      int D5 = x[15];
      int Y6 = x[16];
      int M6 = x[17];
      int D6 = x[18];
      int P = x[19];
      int N = x[20];
      int CpY = x[21];
      int Num1;
      int Den1;
      int Num2;
      int Den2;
      double Frac1;
      double Frac2;
      IntegerVector Num1Diff(6);
      int valuesNum1Diff[] = {Y2,M2,D2,Y3,M3,D3};
      Num1Diff.assign(valuesNum1Diff,valuesNum1Diff+6);
      Num1 = DayDiff(Num1Diff);
      IntegerVector Den1Diff(6);
      int valuesDen1Diff[] = {Y1,M1,D1,Y3,M3,D3};
      Den1Diff.assign(valuesDen1Diff,valuesDen1Diff+6);
      Den1 = DayDiff(Den1Diff);
      IntegerVector Num2Diff(6);
      int valuesNum2Diff[] = {Y4,M4,D4,Y5,M5,D5};
      Num2Diff.assign(valuesNum2Diff,valuesNum2Diff+6);
      Num2 = DayDiff(Num2Diff);
      IntegerVector Den2Diff(6);
      int valuesDen2Diff[] = {Y4,M4,D4,Y6,M6,D6};
      Den2Diff.assign(valuesDen2Diff,valuesDen2Diff+6);
      Den2 = DayDiff(Den2Diff);
      Frac1 = (double)Num1 / (double)Den1;
      Frac2 = (double)Num2 / (double)Den2;
      DIST_value = (Frac1+Frac2+(double)P-(double)N)/(double)CpY;
      // calculating NAccr
      IntegerVector Diff25(6);
      int valuesDiff25[] = {Y2,M2,D2,Y5,M5,D5};
      Diff25.assign(valuesDiff25,valuesDiff25+6);
      NAccr = DayDiff(Diff25);
    } else {
      if (DCC==14) { // Act/365 (Canadian Bond)
        // in this case R should pass a vector of 22 integers:
        // c(DCC,Y1,M1,D1,Y2,M2,D2,Y3,M3,D3,Y4,M4,D4,Y5,M5,D5,Y6,M6,D6,P,N,CpY) with
        // Y1-M1-D1 = PCD(t_a,AD) ; Y2-M2-D2 = t_a ; Y3-M3-D3 = NCD(t_a,AD)
        // Y4-M4-D4 = PCD(t_b,AD) ; Y5-M5-D5 = t_b ; Y6-M6-D6 = NCD(t_b,AD)
        // P = P(t_b,AD) ; N = N(t_a,AD)
        int Y3 = x[7];
        int M3 = x[8];
        int D3 = x[9];
        int Y4 = x[10];
        int M4 = x[11];
        int D4 = x[12];
        int Y5 = x[13];
        int M5 = x[14];
        int D5 = x[15];
        int Y6 = x[16];
        int M6 = x[17];
        int D6 = x[18];
        int P = x[19];
        int N = x[20];
        int CpY = x[21];
        int Num1;
        int Den1;
        int Num2;
        // int Den2;
        double Maple1;
        double Maple2;
        IntegerVector Num1Diff(6);
        int valuesNum1Diff[] = {Y2,M2,D2,Y3,M3,D3};
        Num1Diff.assign(valuesNum1Diff,valuesNum1Diff+6);
        Num1 = DayDiff(Num1Diff);
        IntegerVector Den1Diff(6);
        int valuesDen1Diff[] = {Y1,M1,D1,Y3,M3,D3};
        Den1Diff.assign(valuesDen1Diff,valuesDen1Diff+6);
        Den1 = DayDiff(Den1Diff);
        IntegerVector Num2Diff(6);
        int valuesNum2Diff[] = {Y4,M4,D4,Y5,M5,D5};
        Num2Diff.assign(valuesNum2Diff,valuesNum2Diff+6);
        Num2 = DayDiff(Num2Diff);
        // IntegerVector Den2Diff(6);
        // int valuesDen2Diff[] = {Y4,M4,D4,Y6,M6,D6};
        // Den2Diff.assign(valuesDen2Diff,valuesDen2Diff+6);
        // Den2 = DayDiff(Den2Diff);
        if (Num1==Den1) {
          Maple1 = 1;
        } else {
          if (Num1<((double)365/(double)CpY)) {
            Maple1 = (double)CpY*(double)Num1/(double)365;
          } else {
            int Num3;
            IntegerVector Num3Diff(6);
            int valuesNum3Diff[] = {Y1,M1,D1,Y2,M2,D2};
            Num3Diff.assign(valuesNum3Diff,valuesNum3Diff+6);
            Num3 = DayDiff(Num3Diff);
            Maple1 = (double)1 - (double)CpY*(double)Num3/(double)365;
          }
        }
        if (Num2==0) {
          Maple2 = 0;
        } else {
          if (Num2<(double)365/(double)CpY) {
            Maple2 = (double)CpY*(double)Num2/(double)365;
          } else {
            int Num4;
            IntegerVector Num4Diff(6);
            int valuesNum4Diff[] = {Y5,M5,D5,Y6,M6,D6};
            Num4Diff.assign(valuesNum4Diff,valuesNum4Diff+6);
            Num4 = DayDiff(Num4Diff);
            Maple2 = (double)1 - (double)CpY*(double)Num4/(double)365;
          }
        }
        DIST_value = (Maple1+Maple2+P-N)/CpY;
        // calculating NAccr
        IntegerVector Diff25(6);
        int valuesDiff25[] = {Y2,M2,D2,Y5,M5,D5};
        Diff25.assign(valuesDiff25,valuesDiff25+6);
        NAccr = DayDiff(Diff25);
      } else {
        if (DCC==3) { // Act/Act (AFB)
          if (Y1==Y2) {
            IntegerVector Diff12(6);
            int valuesDiff12[] = {Y1,M1,D1,Y2,M2,D2};
            Diff12.assign(valuesDiff12,valuesDiff12+6);
            NAccr = DayDiff(Diff12);
            if (LeapDayInside(Diff12)==1) {
              DIST_value = (double)NAccr / (double)366;
            } else {
              DIST_value = (double)NAccr / (double)365;
            }
          } else {
            int Date1 = 10000*Y1 + 100*M1 + D1;
            int Date2 = 10000*Y2 + 100*M2 + D2;
            int Date2inY1 = 10000*Y1 + 100*M2 + D2;
            int Date2adj = 10000*(Y1+1) + 100*M2 + D2;
            if ((M2==2)&(D2==29)) {
              if (leap(Y1)==1) {
                Date2inY1 = 10000*Y1 + 100*M2 + D2;
                Date2adj = 10000*(Y1+1) + 100*M2 + 28;
              } else {
                if (leap(Y1+1)==1) {
                  Date2inY1 = 10000*Y1 + 100*M2 + 28;
                  Date2adj = 10000*(Y1+1) + 100*M2 + D2;
                } else {
                  Date2inY1 = 10000*Y1 + 100*M2 + D2;
                  Date2adj = 10000*(Y1+1) + 100*M2 + D2;
                }
              }
            }
            if (Date2adj==Date2) {
              if (Date2inY1<=Date1) {
                IntegerVector Diff12(6);
                int valuesDiff12[] = {Y1,M1,D1,Y2,M2,D2};
                Diff12.assign(valuesDiff12,valuesDiff12+6);
                NAccr = DayDiff(Diff12);
                if (LeapDayInside(Diff12)==1) {
                  DIST_value = (double)NAccr / (double)366;
                } else {
                  DIST_value = (double)NAccr / (double)365;
                }
              } else { // i.e. if (Date1<Date2inY1)
                IntegerVector Diff12(6);
                int valuesDiff12[] = {Y1,M1,D1,Y1,M2,D2};
                Diff12.assign(valuesDiff12,valuesDiff12+6);
                NAccr = DayDiff(Diff12);
                if (LeapDayInside(Diff12)==1) {
                  DIST_value = (double)NAccr / (double)366 + 1;
                } else {
                  DIST_value = (double)NAccr / (double)365 + 1;
                }
                int valuesDiff_for_NAccr[] = {Y1,M2,D2,Y2,M2,D2};
                IntegerVector Diff_for_NAccr(6);
                Diff_for_NAccr.assign(valuesDiff_for_NAccr,valuesDiff_for_NAccr+6);
                int NAccr_Add;
                NAccr_Add = DayDiff(Diff_for_NAccr);
                NAccr = NAccr + NAccr_Add;
              }
            } else { // i.e. if (Date2adj<Date2)
              if (Date2inY1<=Date1) {
                IntegerVector Diff12(6);
                int valuesDiff12[] = {Y1,M1,D1,(Y1+1),M2,D2};
                Diff12.assign(valuesDiff12,valuesDiff12+6);
                NAccr = DayDiff(Diff12);
                if (LeapDayInside(Diff12)==1) {
                  DIST_value = (double)NAccr / (double)366 + (Y2-(Y1+1));
                } else {
                  DIST_value = (double)NAccr / (double)365 + (Y2-(Y1+1));
                }
                int valuesDiff_for_NAccr[] = {(Y1+1),M2,D2,Y2,M2,D2};
                IntegerVector Diff_for_NAccr(6);
                Diff_for_NAccr.assign(valuesDiff_for_NAccr,valuesDiff_for_NAccr+6);
                int NAccr_Add;
                NAccr_Add = DayDiff(Diff_for_NAccr);
                NAccr = NAccr + NAccr_Add;
              } else { // i.e. if (Date1<Date2inY1)
                IntegerVector Diff12(6);
                int valuesDiff12[] = {Y1,M1,D1,Y1,M2,D2};
                Diff12.assign(valuesDiff12,valuesDiff12+6);
                NAccr = DayDiff(Diff12);
                if (LeapDayInside(Diff12)==1) {
                  DIST_value = (double)NAccr / (double)366 + (Y2-Y1);
                } else {
                  DIST_value = (double)NAccr / (double)365 + (Y2-Y1);
                }
                int valuesDiff_for_NAccr[] = {Y1,M2,D2,Y2,M2,D2};
                IntegerVector Diff_for_NAccr(6);
                Diff_for_NAccr.assign(valuesDiff_for_NAccr,valuesDiff_for_NAccr+6);
                int NAccr_Add;
                NAccr_Add = DayDiff(Diff_for_NAccr);
                NAccr = NAccr + NAccr_Add;
              }
            }
          }
        } else {
          if (DCC==4) { // Act/365L; ISMA-Year; Act/365-366
            // in this case R should pass a vector of 9 integers:
            // c(DCC,Y1,M1,D1,Y2,M2,D2,Y3,CpY) with
            // Y1-M1-D1 = t_a ; Y2-M2-D2 = t_b;
            // Y3 = Year figure of the next coupon payment date after t_b
            int Y3 = x[7];
            int CpY = x[8];
            IntegerVector Diff12(6);
            int valuesDiff12[] = {Y1,M1,D1,Y2,M2,D2};
            Diff12.assign(valuesDiff12,valuesDiff12+6);
            NAccr = DayDiff(Diff12);
            if (CpY==1) {
              if (LeapDayInside(Diff12)==1) {
                DIST_value = (double)NAccr / (double)366;
              } else {
                DIST_value = (double)NAccr / (double)365;
              }
            } else { // i.e. CpY != 1
              if (leap(Y3)==1) {
                DIST_value = (double)NAccr / (double)366;
              } else {
                DIST_value = (double)NAccr / (double)365;
              }
            }
          } else {
            if (DCC==15) { // Act/364
              // in this case R should pass a vector of 7 integers:
              // c(DCC,Y1,M1,D1,Y2,M2,D2) with
              // Y1-M1-D1 = t_a ; Y2-M2-D2 = t_b
              IntegerVector Diff12(6);
              int valuesDiff12[] = {Y1,M1,D1,Y2,M2,D2};
              Diff12.assign(valuesDiff12,valuesDiff12+6);
              NAccr = DayDiff(Diff12);
              DIST_value = (double)NAccr / (double)364;
            } else {
              if (DCC==10) { // Act/365; English; DS-Var-BAB = 1
                // in this case R should pass a vector of 7 integers:
                // c(DCC,Y1,M1,D1,Y2,M2,D2) with
                // Y1-M1-D1 = t_a ; Y2-M2-D2 = t_b
                IntegerVector Diff12(6);
                int valuesDiff12[] = {Y1,M1,D1,Y2,M2,D2};
                Diff12.assign(valuesDiff12,valuesDiff12+6);
                NAccr = DayDiff(Diff12);
                DIST_value = (double)NAccr / (double)365;
              } else {
                if (DCC==11) { // Act(NL)/365; Act(No Leap Year)/365
                  // in this case R should pass a vector of 7 integers:
                  // c(DCC,Y1,M1,D1,Y2,M2,D2) with
                  // Y1-M1-D1 = t_a ; Y2-M2-D2 = t_b
                  IntegerVector Diff12(6);
                  int valuesDiff12[] = {Y1,M1,D1,Y2,M2,D2};
                  Diff12.assign(valuesDiff12,valuesDiff12+6);
                  if (LeapDayInside(Diff12)==1) {
                    NAccr = DayDiff(Diff12) - 1;
                  } else {
                    NAccr = DayDiff(Diff12);
                  }
                  DIST_value = (double)NAccr / (double)365;
                } else {
                  if (DCC==12) { // Act/360; French; DS-Var-BAB = 4
                    // in this case R should pass a vector of 7 integers:
                    // c(DCC,Y1,M1,D1,Y2,M2,D2) with
                    // Y1-M1-D1 = t_a ; Y2-M2-D2 = t_b
                    IntegerVector Diff12(6);
                    int valuesDiff12[] = {Y1,M1,D1,Y2,M2,D2};
                    Diff12.assign(valuesDiff12,valuesDiff12+6);
                    NAccr = DayDiff(Diff12);
                    DIST_value = (double)NAccr / (double)360;
                  } else {
                    if (DCC==16) { // Bus/252
                      // in this case R should pass a vector of 8 integers:
                      // c(DCC,Y1,M1,D1,Y2,M2,D2,NonBus) with
                      // Y1-M1-D1 = t_a ; Y2-M2-D2 = t_b and
                      // NonBus = # of non-business days between t_a (incl) and t_b (excl)
                      int NonBus = x[7];
                      IntegerVector Diff12(6);
                      int valuesDiff12[] = {Y1,M1,D1,Y2,M2,D2};
                      Diff12.assign(valuesDiff12,valuesDiff12+6);
                      NAccr = DayDiff(Diff12) - NonBus;
                      DIST_value = (double)NAccr / (double)252;
                    } else {
                      int D1_adj = D1;
                      int D2_adj = D2;
                      if (DCC==5) { // 30/360; Bond Basis; DS-Var-BAB = none
                        // in this case R should pass a vector of 7 integers:
                        // c(DCC,Y1,M1,D1,Y2,M2,D2) with
                        // Y1-M1-D1 = t_a ; Y2-M2-D2 = t_b
                        D1_adj = std::min(D1,30);
                        if ((D2==31)&((D1==30)|(D1==31))) {
                          D2_adj = 30;
                        } else {
                          D2_adj = D2;
                        }
                      } else {
                        if (DCC==6) { // 30E/360; Special German; DS-Var-BAB = 2
                          // in this case R should pass a vector of 7 integers:
                          // c(DCC,Y1,M1,D1,Y2,M2,D2) with
                          // Y1-M1-D1 = t_a ; Y2-M2-D2 = t_b
                          D1_adj = std::min(D1,30);
                          D2_adj = std::min(D2,30);
                        } else {
                          if (DCC==7) { // 30E/360 (ISDA)
                            // in this case R should pass a vector of 10 integers:
                            // c(DCC,Y1,M1,D1,Y2,M2,D2,Y3,M3,D3) with
                            // Y1-M1-D1 = t_a ; Y2-M2-D2 = t_b ; Y3-M3-D3 = t_M
                            int YMat = x[7];
                            int MMat = x[8];
                            int DMat = x[9];
                            int D2isMat;
                            if ((D2==DMat)&(M2==MMat)&(Y2==YMat)) {
                              D2isMat = 1;
                            } else {
                              D2isMat = 0;
                            }
                            IntegerVector Date1(3);
                            IntegerVector Date2(3);
                            int values1[] = {Y1,M1,D1};
                            int values2[] = {Y2,M2,D2};
                            Date1.assign(values1,values1+3);
                            Date2.assign(values2,values2+3);
                            int LDMDate1 = LDM(Date1);
                            int LDMDate2 = LDM(Date2);
                            if (LDMDate1==1) {
                              D1_adj = 30;
                            } else {
                              D1_adj = D1;
                            }
                            if (D2==31) {
                              D2_adj = 30;
                            } else {
                              if ((D2isMat==0)&(M2==2)&(LDMDate2==1)) {
                                D2_adj = 30;
                              } else {
                                D2_adj = D2;
                              }
                            }
                          } else {
                            if (DCC==8) { // 30/360 (German)
                              // in this case R should pass a vector of 7 integers:
                              // c(DCC,Y1,M1,D1,Y2,M2,D2) with
                              // Y1-M1-D1 = t_a ; Y2-M2-D2 = t_b
                              IntegerVector Date1(3);
                              IntegerVector Date2(3);
                              int values1[] = {Y1,M1,D1};
                              int values2[] = {Y2,M2,D2};
                              Date1.assign(values1,values1+3);
                              Date2.assign(values2,values2+3);
                              int LDMDate1 = LDM(Date1);
                              int LDMDate2 = LDM(Date2);
                              if (LDMDate1==1) {
                                D1_adj = 30;
                              } else {
                                D1_adj = D1;
                              }
                              if (LDMDate2==1) {
                                D2_adj = 30;
                              } else {
                                D2_adj = D2;
                              }
                            } else {
                              if ((DCC==9)|      // 30U/360; 30/360 US; DS-Var-BAB = 3
                                  (DCC==13)) {   // 30/365
                                // in this case R should pass a vector of 8 integers:
                                // c(DCC,Y1,M1,D1,Y2,M2,D2,EOM) with
                                // Y1-M1-D1 = t_a ; Y2-M2-D2 = t_b
                                int EOM = x[7];
                                IntegerVector Date1(3);
                                IntegerVector Date2(3);
                                int values1[] = {Y1,M1,D1};
                                int values2[] = {Y2,M2,D2};
                                Date1.assign(values1,values1+3);
                                Date2.assign(values2,values2+3);
                                int D1isLDFeb;
                                int D2isLDFeb;
                                if ((LDM(Date1)==1)&(M1==2)) {
                                  D1isLDFeb = 1;
                                } else {
                                  D1isLDFeb = 0;
                                }
                                if ((LDM(Date2)==1)&(M2==2)) {
                                  D2isLDFeb = 1;
                                } else {
                                  D2isLDFeb = 0;
                                }
                                if ((EOM==1)&(D1isLDFeb==1)) {
                                  D1_adj = 30;
                                } else {
                                  if (D1==31) {
                                    D1_adj = 30;
                                  } else {
                                    D1_adj = D1;
                                  }
                                }
                                if ((EOM==1)&(D1isLDFeb==1)&(D2isLDFeb==1)) {
                                  D2_adj = 30;
                                } else {
                                  if ((D2==31)&((D1==30)|(D1==31))) {
                                    D2_adj = 30;
                                  } else {
                                    D2_adj = D2;
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                      NAccr = 360*(Y2-Y1)+30*(M2-M1)+(D2_adj-D1_adj);
                      if (DCC==13) {
                        DIST_value = (double)NAccr / (double)365;
                      } else {
                        DIST_value = (double)NAccr / (double)360;
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
  double values[] = {NAccr,DIST_value};
  out.assign (values,values+2);   // assigning from array.
  return out;
}


// // This is the outdated DayCountCon-fuction
// // It was replaced by the function DIST in May/June 2018
// // [[Rcpp::export]]
// NumericVector DayCountCon(NumericVector x) {
//   int Y1 = x[0];
//   int M1 = x[1];
//   int D1 = x[2];
//   int Y2 = x[3];
//   int M2 = x[4];
//   int D2 = x[5];
//   int Y3 = x[6];
//   int M3 = x[7];
//   int D3 = x[8];
//   int CpY = x[9];
//   int DCC = x[10];
//   double NAccr;
//   double NPeriod;
//   double AccrFactor;
//   double CoupFactor;
//   NumericVector out(4);
//   if (DCC==1) { // Act/Act ISDA; DS-Var-BAB = 6
//     IntegerVector Diff12(6);
//     IntegerVector Diff13(6);
//     int valuesDiff12[] = {Y1,M1,D1,Y2,M2,D2};
//     int valuesDiff13[] = {Y1,M1,D1,Y3,M3,D3};
//     Diff12.assign(valuesDiff12,valuesDiff12+6);
//     Diff13.assign(valuesDiff13,valuesDiff13+6);
//     int Ydiff21 = Y2 - Y1;
//     int Ydiff31 = Y3 - Y1;
//     int DaysInY1 = DaysInYear(Y1);
//     if (Ydiff21==0) {
//       NAccr = DayDiff(Diff12);
//       AccrFactor = (double)NAccr / (double)DaysInY1;
//     } else {
//       int lengthYVect21 = Ydiff21 + 1;
//       IntegerVector YVector(lengthYVect21);
//       for(int i = 0; i < lengthYVect21; ++i) {
//         YVector[i] = DaysInYear(Y1 + i);
//       }
//       int valuesDaysY1[] = {Y1,M1,D1,Y1+1,1,1};
//       IntegerVector DaysY1Vector(6);
//       DaysY1Vector.assign(valuesDaysY1,valuesDaysY1+6);
//       int DaysY1 = DayDiff(DaysY1Vector);
//       int valuesDaysY2[] = {Y2,1,1,Y2,M2,D2};
//       IntegerVector DaysY2Vector(6);
//       DaysY2Vector.assign(valuesDaysY2,valuesDaysY2+6);
//       int DaysY2 = DayDiff(DaysY2Vector);
//       IntegerVector AccrDayVec(lengthYVect21);
//       AccrDayVec[0] = DaysY1;
//       for(int i = 1; i < lengthYVect21-1; ++i) {
//         AccrDayVec[i] = YVector[i];
//       }
//       AccrDayVec[lengthYVect21-1] = DaysY2;
//       NAccr = accumulate(AccrDayVec.begin(), AccrDayVec.end(), 0);
//       NumericVector AccrFacVec(lengthYVect21);
//       for (int i = 0; i < lengthYVect21; ++i) {
//         double Num = AccrDayVec[i];
//         double Den = YVector[i];
//         AccrFacVec[i] = Num / Den;
//       }
//       AccrFactor = sumC(AccrFacVec);
//     }
//     if (Ydiff31==0) {
//       NPeriod = DayDiff(Diff13);
//       CoupFactor = (double)NPeriod / (double)DaysInY1;
//     } else {
//       int lengthYVect31 = Ydiff31 + 1;
//       IntegerVector YVector(lengthYVect31);
//       for(int i = 0; i < lengthYVect31; ++i) {
//         YVector[i] = DaysInYear(Y1 + i);
//       }
//       int valuesDaysY1[] = {Y1,M1,D1,Y1+1,1,1};
//       IntegerVector DaysY1Vector(6);
//       DaysY1Vector.assign(valuesDaysY1,valuesDaysY1+6);
//       int DaysY1 = DayDiff(DaysY1Vector);
//       int valuesDaysY3[] = {Y3,1,1,Y3,M3,D3};
//       IntegerVector DaysY3Vector(6);
//       DaysY3Vector.assign(valuesDaysY3,valuesDaysY3+6);
//       int DaysY3 = DayDiff(DaysY3Vector);
//       IntegerVector PeriodDayVec(lengthYVect31);
//       PeriodDayVec[0] = DaysY1;
//       for(int i = 1; i < lengthYVect31-1; ++i) {
//         PeriodDayVec[i] = YVector[i];
//       }
//       PeriodDayVec[lengthYVect31-1] = DaysY3;
//       NPeriod = accumulate(PeriodDayVec.begin(), PeriodDayVec.end(), 0);
//       NumericVector PeriodFacVec(lengthYVect31);
//       for (int i = 0; i < lengthYVect31; ++i) {
//         double Num = PeriodDayVec[i];
//         double Den = YVector[i];
//         PeriodFacVec[i] = Num / Den;
//       }
//       CoupFactor = sumC(PeriodFacVec);
//     }
//   } else {
//     if (DCC==2) { // Act/Act (ICMA); Act/Act (ISMA); DS-Var-BAB = 5
//       IntegerVector Diff12(6);
//       IntegerVector Diff13(6);
//       int valuesDiff12[] = {Y1,M1,D1,Y2,M2,D2};
//       int valuesDiff13[] = {Y1,M1,D1,Y3,M3,D3};
//       Diff12.assign(valuesDiff12,valuesDiff12+6);
//       Diff13.assign(valuesDiff13,valuesDiff13+6);
//       NAccr = DayDiff(Diff12);
//       NPeriod = DayDiff(Diff13);
//       AccrFactor = (double)NAccr / ((double)NPeriod * (double)CpY);
//       CoupFactor = (double)1 / (double)CpY;
//     } else {
//       if (DCC==3) { // Act/365; English; DS-Var-BAB = 1
//         IntegerVector Diff12(6);
//         IntegerVector Diff13(6);
//         int valuesDiff12[] = {Y1,M1,D1,Y2,M2,D2};
//         int valuesDiff13[] = {Y1,M1,D1,Y3,M3,D3};
//         Diff12.assign(valuesDiff12,valuesDiff12+6);
//         Diff13.assign(valuesDiff13,valuesDiff13+6);
//         NAccr = DayDiff(Diff12);
//         NPeriod = DayDiff(Diff13);
//         AccrFactor = (double)NAccr / (double)365;
//         CoupFactor = (double)NPeriod/ (double)365;
//       } else {
//         if (DCC==4) { // Act/360; French; DS-Var-BAB = 4
//           IntegerVector Diff12(6);
//           IntegerVector Diff13(6);
//           int valuesDiff12[] = {Y1,M1,D1,Y2,M2,D2};
//           int valuesDiff13[] = {Y1,M1,D1,Y3,M3,D3};
//           Diff12.assign(valuesDiff12,valuesDiff12+6);
//           Diff13.assign(valuesDiff13,valuesDiff13+6);
//           NAccr = DayDiff(Diff12);
//           NPeriod = DayDiff(Diff13);
//           AccrFactor = (double)NAccr / (double)360;
//           CoupFactor = (double)NPeriod / (double)360;
//         } else {
//           if (DCC==5) { // 30/360; Bond Basis; DS-Var-BAB = none
//             D1 = std::min(D1,30);
//             if (D1==30.0) {
//               D2 = std::min(D2,30);
//               D3 = std::min(D3,30);
//             }
//           } else {
//             if (DCC==6) { // 30E/360; Special German; DS-Var-BAB = 2
//               D1 = std::min(D1,30);
//               D2 = std::min(D2,30);
//               D3 = std::min(D3,30);
//             } else {
//               if (DCC==7) { // 30E/360 ISDA; German; DS-Var-BAB = none
//                 int YMat = x[11];
//                 int MMat = x[12];
//                 int DMat = x[13];
//                 IntegerVector Date1(3);
//                 IntegerVector Date2(3);
//                 IntegerVector Date3(3);
//                 int values1[] = {Y1,M1,D1};
//                 int values2[] = {Y2,M2,D2};
//                 int values3[] = {Y3,M3,D3};
//                 Date1.assign(values1,values1+3);
//                 Date2.assign(values2,values2+3);
//                 Date3.assign(values3,values3+3);
//                 int LDMDate1 = LDM(Date1);
//                 int LDMDate2 = LDM(Date2);
//                 int LDMDate3 = LDM(Date3);
//                 if (LDMDate1==1) {
//                   D1 = 30;
//                 }
//                 if (LDMDate2==1) {
//                   if (M2==2) {
//                     if (!((D2==DMat)&(M2==MMat)&(Y2==YMat))) {
//                       D2 = 30;
//                     }
//                   } else {
//                     D2 = 30;
//                   }
//                 }
//                 if (LDMDate3==1) {
//                   if (M3==2) {
//                     if (!((D3==DMat)&(M3==MMat)&(Y3==YMat))) {
//                       D3 = 30;
//                     }
//                   } else {
//                     D3 = 30;
//                   }
//                 }
//               } else {
//                 if (DCC==8) { // 30U/360; 30/360 US; DS-Var-BAB = 3
//                   IntegerVector Date1(3);
//                   IntegerVector Date2(3);
//                   IntegerVector Date3(3);
//                   int values1[] = {Y1,M1,D1};
//                   int values2[] = {Y2,M2,D2};
//                   int values3[] = {Y3,M3,D3};
//                   Date1.assign(values1,values1+3);
//                   Date2.assign(values2,values2+3);
//                   Date3.assign(values3,values3+3);
//                   int LDMDate1 = LDM(Date1);
//                   int LDMDate2 = LDM(Date2);
//                   int LDMDate3 = LDM(Date3);
//                   if (((M2==2)&(LDMDate2==1))&((M1==2)&(LDMDate1=1))) {
//                     D2 = 30;
//                   }
//                   if (((M3==2)&(LDMDate3==1))&((M1==2)&(LDMDate1=1))) {
//                     D3 = 30;
//                   }
//                   if ((M1==2)&(LDMDate1=1)) {
//                     D1 = 30;
//                   }
//                   if ((D2==31)&((D1==30)|(D1==31))) {
//                     D2 = 30;
//                   }
//                   if ((D3==31)&((D1==30)|(D1==31))) {
//                     D3 = 30;
//                   }
//                   if (D1==31) {
//                     D1 = 30;
//                   }
//                 }
//               }
//             }
//           }
//           NAccr = 360*(Y2-Y1)+30*(M2-M1)+(D2-D1);
//           NPeriod = 360*(Y3-Y1)+30*(M3-M1)+(D3-D1);
//           AccrFactor = (double)NAccr / (double)360;
//           CoupFactor = (double)NPeriod / (double)360;
//         }
//       }
//     }
//   }
//   double values[] = {NAccr,NPeriod,AccrFactor,CoupFactor};
//   out.assign (values,values+4);   // assigning from array.
//   return out;
// }



// [[Rcpp::export]]
NumericVector PayCalc(NumericMatrix x) {
  // this function uses the function DIST
  //
  // for DCC = {1,3,5,6,8,10,11} x is a NumericMatrix with 9 columns,
  // each row containing the vector:
  // c(RV,Coup,DCC,Y1,M1,D1,Y2,M2,D2) with
  // Y1-M1-D1 = t_a ; Y2-M2-D2 = t_b
  //
  // for DCC = 2 x is a NumericMatrix with 24 columns,
  // each row containing the vector:
  // c(RV,Coup,DCC,Y1,M1,D1,Y2,M2,D2,Y3,M3,D3,Y4,M4,D4,Y5,M5,D5,Y6,M6,D6,P,N,CpY) with
  // Y1-M1-D1 = PCD(t_a,AD) ; Y2-M2-D2 = t_a ; Y3-M3-D3 = NCD(t_a,AD)
  // Y4-M4-D4 = PCD(t_b,AD) ; Y5-M5-D5 = t_b ; Y6-M6-D6 = NCD(t_b,AD)
  // P = P(t_b,AD) ; N = N(t_a,AD)
  //
  // for DCC = 4 x is a NumericMatrix of 11 columns,
  // each row containing the vector:
  // c(RV,Coup,DCC,Y1,M1,D1,Y2,M2,D2,Y3,CpY) with
  // Y1-M1-D1 = t_a ; Y2-M2-D2 = t_b; Y3-M3-D3 = t_c
  //
  // for DCC = 7 x is a NumericMatrix with 12 columns,
  // each row containing the vector:
  // c(RV,Coup,DCC,Y1,M1,D1,Y2,M2,D2,Y3,M3,D3) with
  // Y1-M1-D1 = t_a ; Y2-M2-D2 = t_b ; Y3-M3-D3 = t_M
  //
  // for DCC = 9 x is a NumericMatrix with 10 columns,
  // each row containing the vector:
  // c(RV,Coup,DCC,Y1,M1,D1,Y2,M2,D2,EOM) with
  // Y1-M1-D1 = t_a ; Y2-M2-D2 = t_b
  //
  int nrow = x.nrow();
  int ncol = x.ncol();
  NumericVector out(nrow);
  for (int i = 0; i < nrow; i++) {
    double RV = x(i,0);
    double Coup = x(i,1);
    NumericMatrix DIST_Input = x(Range(i,i),Range(2,(ncol-1)));
    double CoupFactor = DIST(DIST_Input)[1];
    double CoupPayment = RV * Coup * CoupFactor;
    out[i] = CoupPayment;
  }
  return out;
}


// This is the outdated PayCalc-fuction
// It was replaced by a new version in May/June 2018
// // [[Rcpp::export]]
// NumericVector PayCalc(NumericMatrix x) {
//   int nrow = x.nrow();
//   NumericVector out(nrow);
//   for (int i = 0; i < nrow; i++) {
//     int Y1 = x(i,0);
//     int M1 = x(i,1);
//     int D1 = x(i,2);
//     int Y2 = x(i,3);
//     int M2 = x(i,4);
//     int D2 = x(i,5);
//     int Y3 = x(i,6);
//     int M3 = x(i,7);
//     int D3 = x(i,8);
//     int CpY = x(i,9);
//     int DCC = x(i,10);
//     int YMat = x(i,11);
//     int MMat = x(i,12);
//     int DMat = x(i,13);
//     double RV = x(i,14);
//     double Coup = x(i,15);
//     NumericVector Input(14);
//     int InputValues[] = {Y1,M1,D1,Y2,M2,D2,Y3,M3,D3,CpY,DCC,YMat,MMat,DMat};
//     Input.assign (InputValues,InputValues+14);
//     NumericVector DayCountOut = DayCountCon(Input);
//     double CoupFactor = DayCountOut[3];
//     double CoupPayment = RV * Coup * CoupFactor;
//     out[i] = CoupPayment;
//   }
//   return out;
// }


// // This is the outdated AccrIntCalc-fuction
// // Its functionality was replaced by DIST in July 2018
// // [[Rcpp::export]]
// NumericVector AccrIntCalc(NumericMatrix x) {
//   int nrow = x.nrow();
//   NumericVector out(nrow);
//   for (int i = 0; i < nrow; i++) {
//     int Y1 = x(i,0);
//     int M1 = x(i,1);
//     int D1 = x(i,2);
//     int Y2 = x(i,3);
//     int M2 = x(i,4);
//     int D2 = x(i,5);
//     int Y3 = x(i,6);
//     int M3 = x(i,7);
//     int D3 = x(i,8);
//     int CpY = x(i,9);
//     int DCC = x(i,10);
//     int YMat = x(i,11);
//     int MMat = x(i,12);
//     int DMat = x(i,13);
//     double RV = x(i,14);
//     double Coup = x(i,15);
//     NumericVector Input(14);
//     int InputValues[] = {Y1,M1,D1,Y2,M2,D2,Y3,M3,D3,CpY,DCC,YMat,MMat,DMat};
//     Input.assign (InputValues,InputValues+14);
//     NumericVector DayCountOut = DayCountCon(Input);
//     double AccrFactor = DayCountOut[2];
//     double AccrInt = RV * Coup * AccrFactor;
//     out[i] = AccrInt;
//   }
//   return out;
// }



// [[Rcpp::export]]
IntegerVector NumToDate(IntegerVector x) {
  // Only the first 4 elements of x are used.
  // The first element is assumed to be the number of days from origin that has to be converted to a date.
  // Elements 2 to 4 are optional and can represent the origin of the first element as a vector (yyyy,mm,dd).
  // If x consists of 1 element only, the origin is by default 1970-01-01.
  int INPUT = x[0];
  if (x.size()==4) {
    int OrigYear = x[1];
    int OrigMonth = x[2];
    int OrigDay = x[3];
    IntegerVector IN_adj;
    if (OrigYear>1969) {
      int IN_adj_values[] = {1970,1,1,OrigYear,OrigMonth,OrigDay};
      IN_adj.assign(IN_adj_values,IN_adj_values+6);
      INPUT = INPUT + DayDiff(IN_adj);
    } else {
      int IN_adj_values[] = {OrigYear,OrigMonth,OrigDay,1970,1,1};
      IN_adj.assign(IN_adj_values,IN_adj_values+6);
      INPUT = INPUT - DayDiff(IN_adj);
    }
  }
  int YEAR;
  int MONTH = 0;
  int DAY = 0;
  int Rest;
  IntegerVector MonthDays(13);
  int InputValues[] = {0,31,59,90,120,151,181,212,243,273,304,334,365};      // accumulated days in non-leap year
  MonthDays.assign(InputValues,InputValues+13);
  IntegerVector MonthDaysLeap(13);
  int InputValuesLeap[] = {0,31,60,91,121,152,182,213,244,274,305,335,366};  // accumulated days in leap year
  MonthDaysLeap.assign(InputValuesLeap,InputValuesLeap+13);
  if (INPUT==0) {
    YEAR = 1970;
    MONTH = 1;
    DAY = 1;
  } else {
    if (INPUT>0) {
      if (INPUT<365) {
        YEAR = 1970;
        while (MonthDays[MONTH] < (INPUT+1)) {
          MONTH++;                                                           // finding the month
        }
        DAY = (INPUT+1) - MonthDays[MONTH - 1];                              // finding the day
      } else {
        if (INPUT<730) {
          YEAR = 1971;
          INPUT = INPUT - 365;
          while (MonthDays[MONTH] < (INPUT+1)) {
            MONTH++;                                                         // finding the month
          }
          DAY = (INPUT+1) - MonthDays[MONTH - 1];                            // finding the day
        } else {
          // here vector stuff can start since vector length will be at least 2
          int ApproxNumYear = round((double)(INPUT)/(double)365.25);            // this is the approximate number of years covered by the input INPUT
          IntegerVector HelpYear = seq_len(ApproxNumYear);                      // creating the vector HelpYear = (1,2,3,4,...,ApproxNumYear)
          IntegerVector Years(ApproxNumYear);                                   // creating the vector Years of length ApproxNumYear
          Years = 1969 + HelpYear;                                              // modifying Years to Years = (1970,1971,1972,...,(1969+ApproxNumYear))
          IntegerVector DaysInYears(ApproxNumYear);                             // creating the vector DaysInYears of length ApproxNumYear
          for(int i = 0; i < ApproxNumYear; ++i) {                              // filling DaysInYears with the respective year's number of days
            DaysInYears[i] = DaysInYear(Years[i]);                              // DaysInYears = (365,365,366,365,...,DaysInYear(1969+ApproxNumYear))
          }
          int Ysum = accumulate(DaysInYears.begin(), DaysInYears.end(), 0);     // summing up the elements of DaysInYears
          Rest = INPUT - Ysum;
          if (Rest==0) {
            YEAR = Years[ApproxNumYear-1]+1;
            MONTH = 1;
            DAY = 1;
          } else {
            if (Rest<0) {
              YEAR = Years[ApproxNumYear-1];
              DaysInYears[ApproxNumYear-1] = 0;
              Ysum = accumulate(DaysInYears.begin(), DaysInYears.end(), 0);
              Rest = INPUT - Ysum;
            } else {
              YEAR = Years[ApproxNumYear-1]+1;
            }
            int LYLeap = leap(YEAR);
            if (LYLeap==1) {
              while (MonthDaysLeap[MONTH] < (Rest+1)) {
                MONTH++;                                                               // finding the month
              }
              DAY = (Rest+1) - MonthDaysLeap[MONTH - 1];                                   // finding the day
            } else {
              while (MonthDays[MONTH] < (Rest+1)) {
                MONTH++;                                                               // finding the month
              }
              DAY = (Rest+1) - MonthDays[MONTH - 1];                                   // finding the day
            }
          }
        }
      }
    } else {
      int absINPUT = (-1)*INPUT;
      if (absINPUT<366) {
        YEAR = 1969;
        absINPUT = 365 - absINPUT;
        if (absINPUT==0) {
          MONTH = 1;
          DAY = 1;
        } else {
          while (MonthDays[MONTH] < (absINPUT+1)) {
            MONTH++;                                                              // finding the month
          }
          DAY = (absINPUT+1) - MonthDays[MONTH - 1];                              // finding the day
        }
      } else {
        if (absINPUT<731) {
          YEAR = 1968;
          absINPUT = 731 - absINPUT;
          while (MonthDaysLeap[MONTH] < (absINPUT+1)) {
            MONTH++;                                                              // finding the month
          }
          DAY = (absINPUT+1) - MonthDaysLeap[MONTH - 1];                              // finding the day
        } else {
          // here vector stuff can start since vector length will be at least 2
          int ApproxNumYear = round((double)(absINPUT)/(double)365.25);                // this is the approximate number of years covered by the input absabsINPUT
          IntegerVector HelpYear = seq_len(ApproxNumYear);          // creating the vector HelpYear = (1,2,3,4,...,ApproxNumYear)
          IntegerVector Years(ApproxNumYear);                       // creating the vector Years of length ApproxNumYear
          Years = 1970 - HelpYear;                                  // modifying Years to Years = (1970,1971,1972,...,(1969+ApproxNumYear))
          IntegerVector DaysInYears(ApproxNumYear);                 // creating the vector DaysInYears of length ApproxNumYear
          for(int i = 0; i < ApproxNumYear; ++i) {                  // filling DaysInYears with the respective year's number of days
            DaysInYears[i] = DaysInYear(Years[i]);                  // DaysInYears = (365,365,366,365,...,DaysInYear(1969+ApproxNumYear))
          }
          int Ysum = accumulate(DaysInYears.begin(), DaysInYears.end(), 0);       // summing up the elements of DaysInYears
          Rest = Ysum - absINPUT;
          if (Rest==0) {
            YEAR = Years[ApproxNumYear-1];                          // the year that corresponds to absINPUT
            MONTH = 1;
            DAY = 1;
          } else {
            if (Rest<0) {
              YEAR = Years[ApproxNumYear-1] - 1;
              Ysum = Ysum + DaysInYear(YEAR);
              Rest = Ysum - absINPUT;
            } else {
              YEAR = Years[ApproxNumYear-1];
            }
            int LYLeap = leap(YEAR);
            if (LYLeap==1) {
              while (MonthDaysLeap[MONTH] < (Rest+1)) {
                MONTH++;                                                               // finding the month
              }
              DAY = (Rest+1) - MonthDaysLeap[MONTH - 1];                                   // finding the day
            } else {
              while (MonthDays[MONTH] < (Rest+1)) {
                MONTH++;                                                               // finding the month
              }
              DAY = (Rest+1) - MonthDays[MONTH - 1];                                   // finding the day
            }
          }
        }
      }
    }
  }
  IntegerVector out;
  int outValues[] = {YEAR,MONTH,DAY};
  out.assign(outValues,outValues+3);
  return out;
}



// [[Rcpp::export]]
IntegerVector CppPrevDate(IntegerVector x) {
  int Atom1AD1 = x[0];
  int Atom2AD1 = x[1];
  int Atom3AD1 = x[2];
  int Atom1Em = x[3];
  int Atom2Em = x[4];
  int Atom3Em = x[5];
  int Atom1Refer = x[6];
  int Atom2Refer = x[7];
  int Atom3Refer = x[8];
  int CpY = x[9];
  int EOM = x[10];
  IntegerVector out(3);
  int gap = round(365/CpY);
  IntegerVector Diff_AD1_Em(6);
    int Diff_AD1_Em_values[] = {Atom1AD1,Atom2AD1,Atom3AD1,Atom1Em,Atom2Em,Atom3Em};
    Diff_AD1_Em.assign(Diff_AD1_Em_values,Diff_AD1_Em_values+6);
  int AD1_to_Em = DayDiff(Diff_AD1_Em);
  if (!(AD1_to_Em>0)) {
    IntegerVector outConstruction(6);
      int outConstruction_values[] = {1970,1,1,Atom1AD1,Atom2AD1,15};
      outConstruction.assign(outConstruction_values,outConstruction_values+6);
    int out_Num = DayDiff(outConstruction) - gap;
    IntegerVector out_conv(4);
      int out_conv_values[] = {out_Num,1970,1,1};
      out_conv.assign(out_conv_values,out_conv_values+4);
    out = NumToDate(out_conv);
    IntegerVector Refer(3);
      int Refer_values[] = {Atom1Refer,Atom2Refer,Atom3Refer};
      Refer.assign(Refer_values,Refer_values+3);
    // int LDM_Refer = LDM(Refer);
    // if (LDM_Refer==1) {
    if (EOM==1) {
      out(2) = DaysInMonth(out);
    } else {
      IntegerVector out_A(3);
        int out_A_values[] = {out(0),out(1),Atom3Refer};
        out_A.assign(out_A_values,out_A_values+3);
      if (out_A(2)>DaysInMonth(out_A)) {
        out(2) = DaysInMonth(out);
      } else {
        out = out_A;
      }
    }
  } else {
    int out_values[] = {NA_INTEGER,NA_INTEGER,NA_INTEGER};
    out.assign(out_values,out_values+3);
  }
  return out;
}



// [[Rcpp::export]]
IntegerVector CppSuccDate(IntegerVector x) {
  int Atom1ADfin = x[0];
  int Atom2ADfin = x[1];
  int Atom3ADfin = x[2];
  int Atom1Mat = x[3];
  int Atom2Mat = x[4];
  int Atom3Mat = x[5];
  int Atom1Refer = x[6];
  int Atom2Refer = x[7];
  int Atom3Refer = x[8];
  int CpY = x[9];
  int EOM = x[10];
  IntegerVector out(3);
  int gap = round(365/CpY);
  IntegerVector Diff_Mat_ADfin(6);
  int Diff_Mat_ADfin_values[] = {Atom1Mat,Atom2Mat,Atom3Mat,Atom1ADfin,Atom2ADfin,Atom3ADfin};
  Diff_Mat_ADfin.assign(Diff_Mat_ADfin_values,Diff_Mat_ADfin_values+6);
  int Mat_to_ADfin = DayDiff(Diff_Mat_ADfin);
  if (!(Mat_to_ADfin>0)) {
    IntegerVector outConstruction(6);
    int outConstruction_values[] = {1970,1,1,Atom1ADfin,Atom2ADfin,15};
    outConstruction.assign(outConstruction_values,outConstruction_values+6);
    int out_Num = DayDiff(outConstruction) + gap;
    IntegerVector out_conv(4);
    int out_conv_values[] = {out_Num,1970,1,1};
    out_conv.assign(out_conv_values,out_conv_values+4);
    out = NumToDate(out_conv);
    IntegerVector Refer(3);
    int Refer_values[] = {Atom1Refer,Atom2Refer,Atom3Refer};
    Refer.assign(Refer_values,Refer_values+3);
    // int LDM_Refer = LDM(Refer);
    // if (LDM_Refer==1) {
    if (EOM==1) {
      out(2) = DaysInMonth(out);
    } else {
      IntegerVector out_A(3);
      int out_A_values[] = {out(0),out(1),Atom3Refer};
      out_A.assign(out_A_values,out_A_values+3);
      if (out_A(2)>DaysInMonth(out_A)) {
        out(2) = DaysInMonth(out);
      } else {
        out = out_A;
      }
    }
  } else {
    int out_values[] = {NA_INTEGER,NA_INTEGER,NA_INTEGER};
    out.assign(out_values,out_values+3);
  }
  return out;
}






// // [[Rcpp::export]]
// NumericVector callFunction(NumericVector x, Function f) {
//   NumericVector res = f(x);
//   return res;
// }

// // [[Rcpp::export]]
// NumericVector NewtonRaphson(Function func, Function d_func, double Start, double Prec) {
//   NumericVector out;
//   int i = 0;
//   double v = Start;
//   double f_value;
//   double d_fvalue;
//   f_value = as<double>(func(v));
//   d_fvalue = as<double>(d_func(v));
//   double vNext = v - f_value/d_fvalue;
//   while(abs(vNext-v)>Prec) {
//     ++i;
//     v = vNext;
//     f_value = as<double>(func(v));
//     d_fvalue = as<double>(d_func(v));
//     vNext = v - f_value/d_fvalue;
//   }
//   double values[] = {vNext,(double)i,f_value};
//   out.assign (values,values+3);   // assigning from array.
//   return out;
// }



/*** R


#_______________________________________________________________________________________________________________________________________
#***************************************************************************************************************************************
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
########################################################################################################################################
#######------>>>                                       <<<------########################################################################
####------->>>   testing and benchmarking CppSuccDate    <<<------######################################################################
#######------>>>                                       <<<------########################################################################
#
#
#
#################################################################################################################################
######                                                                    #######################################################
######   testing whether RSuccDate and CppSuccDate do exactly the same    #######################################################
#
# library(timeDate)
#
# years1<-sample(c(1000:3000),1000,replace = TRUE)
# years2<-sample(c(1000:3000),1000,replace = TRUE)
# years3<-sample(c(1000:3000),1000,replace = TRUE)
# months1<-sample(c(1:12),1000,replace = TRUE)
# months2<-sample(c(1:12),1000,replace = TRUE)
# months3<-sample(c(1:12),1000,replace = TRUE)
# days1<-sample(c(1:28),1000,replace = TRUE)
# days2<-sample(c(1:28),1000,replace = TRUE)
# days3<-sample(c(1:28),1000,replace = TRUE)
# CpY<-sample(c(1,2,3,4,6,12),1000,replace = TRUE)
# SuccDateTest<-data.frame(years1,months1,days1,years2,months2,days2,years3,months3,days3,CpY)
# SuccDateTest$Date1<-as.Date(gsub("[ //s]","",paste(SuccDateTest$years1,"-",SuccDateTest$months1,"-",SuccDateTest$days1)))
# SuccDateTest$Date2<-as.Date(gsub("[ //s]","",paste(SuccDateTest$years2,"-",SuccDateTest$months2,"-",SuccDateTest$days2)))
# SuccDateTest$Date3<-as.Date(gsub("[ //s]","",paste(SuccDateTest$years3,"-",SuccDateTest$months3,"-",SuccDateTest$days3)))
#
# # creating the anniversary date preceding Em in R
# RSuccDate<-function(AnnivDatesfin,Mat,Refer,CpY) {
#   gap<-round(abs(as.numeric(difftime(as.Date(ISOdatetime(t(atoms(as.timeDate(Mat))[1]),t(atoms(as.timeDate(Mat))[2]),15,12,0,0)),
#                                      as.Date(ISOdatetime(t(atoms(as.timeDate(Mat))[1]+1),t(atoms(as.timeDate(Mat))[2]),15,12,0,0)))))/CpY,digits=0)
#   if (!(Mat<AnnivDatesfin)) {
#     SuccDate<-as.Date(ISOdatetime(atoms(as.timeDate(AnnivDatesfin))[1],atoms(as.timeDate(AnnivDatesfin))[2],15,12,0,0))
#     SuccDate<-SuccDate+gap
#     if ((as.Date(timeLastDayInMonth(Refer))==Refer)==TRUE) {
#       SuccDate<-as.Date(timeLastDayInMonth(SuccDate))
#     } else {
#       SuccDate_A<-as.Date(ISOdatetime(t(atoms(as.timeDate(SuccDate))[1]),t(atoms(as.timeDate(SuccDate))[2]),t(atoms(as.timeDate(Refer))[3]),12,0,0))
#       if (is.na(SuccDate_A)) {
#         SuccDate<-as.Date(timeLastDayInMonth(SuccDate))
#       } else {
#         SuccDate<-SuccDate_A
#       }
#     }
#   } else {
#     SuccDate<-as.Date(NA)
#   }
#   SuccDate<-as.Date(SuccDate)
#   return(SuccDate)
# }
#
# # creating the anniversary date preceding Em in C++
# CSuccDate<-function(x) {
#   out<-CppSuccDate(x)
#   if (any(is.na(out))) {
#     out<-as.Date(NA)
#   } else {
#     out<-as.Date(paste(out,collapse="-"),origin="1970-01-01")
#   }
#   return(out)
# }
#
# for (i in c(1:nrow(SuccDateTest))) {
#   SuccDateTest$RSuccDate[i]<-RSuccDate(SuccDateTest$Date1[i],SuccDateTest$Date2[i],SuccDateTest$Date3[i],SuccDateTest$CpY[i])
# }
#
# for (i in c(1:nrow(SuccDateTest))) {
#   SuccDateTest$CSuccDate[i]<-CSuccDate(c(SuccDateTest$years1[i],SuccDateTest$months1[i],SuccDateTest$days1[i],
#                                 SuccDateTest$years2[i],SuccDateTest$months2[i],SuccDateTest$days2[i],
#                                 SuccDateTest$years3[i],SuccDateTest$months3[i],SuccDateTest$days3[i],
#                                 SuccDateTest$CpY[i]))
# }
#
# for (i in c(1:nrow(SuccDateTest))) {
#   if (all(is.na(SuccDateTest$CSuccDate[i]),is.na(SuccDateTest$RSuccDate[i]))) {
#     SuccDateTest$Check[i]<-1
#   } else {
#     if (SuccDateTest$CSuccDate[i]==SuccDateTest$RSuccDate[i]) {
#       SuccDateTest$Check[i]<-1
#     } else {
#       SuccDateTest$Check[i]<-0
#     }
#   }
# }
#
# length(which(SuccDateTest$Check==0))
#
#
######   HURRRAAAAYYYY!!! RSuccDate and CSuccDate do exactly the same!    #######################################################
######                                                                    #######################################################
#################################################################################################################################
###^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^###
###vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv###
#################################################################################################################################
######                                            ###############################################################################
######   benchmarking RSuccDate and CSuccDate     ###############################################################################
#
# years1<-sample(c(1000:3000),100,replace = TRUE)
# years2<-sample(c(1000:3000),100,replace = TRUE)
# years3<-sample(c(1000:3000),100,replace = TRUE)
# months1<-sample(c(1:12),100,replace = TRUE)
# months2<-sample(c(1:12),100,replace = TRUE)
# months3<-sample(c(1:12),100,replace = TRUE)
# days1<-sample(c(1:28),100,replace = TRUE)
# days2<-sample(c(1:28),100,replace = TRUE)
# days3<-sample(c(1:28),100,replace = TRUE)
# CpY<-sample(c(1,2,3,4,6,12),100,replace = TRUE)
# SuccDateTest<-data.frame(years1,months1,days1,years2,months2,days2,years3,months3,days3,CpY)
# SuccDateTest$Date1<-as.Date(gsub("[ //s]","",paste(SuccDateTest$years1,"-",SuccDateTest$months1,"-",SuccDateTest$days1)))
# SuccDateTest$Date2<-as.Date(gsub("[ //s]","",paste(SuccDateTest$years2,"-",SuccDateTest$months2,"-",SuccDateTest$days2)))
# SuccDateTest$Date3<-as.Date(gsub("[ //s]","",paste(SuccDateTest$years3,"-",SuccDateTest$months3,"-",SuccDateTest$days3)))
#
#
# library(rbenchmark)
#
# DF_RSuccDate<-function(x) {
#   OUT<-RSuccDate(x[1,1],x[1,2],x[1,3],x[1,4])
#   for (i in c(2:nrow(x))) {
#     OUTnext<-RSuccDate(x[i,1],x[i,2],x[i,3],x[i,4])
#     OUT<-append(OUT,OUTnext)
#   }
#   return(OUT)
# }
#
# DF_CSuccDate<-function(x) {
#   OUT<-CSuccDate(c(x[1,1],x[1,2],x[1,3],x[1,4],x[1,5],x[1,6],x[1,7],x[1,8],x[1,9],x[1,10]))
#   for (i in c(2:nrow(x))) {
#     OUTnext<-CSuccDate(c(x[i,1],x[i,2],x[i,3],x[i,4],x[i,5],x[i,6],x[i,7],x[i,8],x[i,9],x[i,10]))
#     OUT<-append(OUT,OUTnext)
#   }
#   return(OUT)
# }
#
# x<-SuccDateTest[,c(11:13,10)]
# y<-SuccDateTest[,c(1:10)]
#
# # check that DF_RSuccDate is the same as DF_CSuccDate
# stopifnot(all.equal(DF_RSuccDate(x), DF_CSuccDate(y)))
#
# benchmark(DF_RSuccDate(x), DF_CSuccDate(y), order="relative")[,1:4]
#
# # benchmarking CSuccDate and RSuccDate with a data frame with 1000 rows
# #              test replications elapsed relative
# # 2 DF_CSuccDate(y)          100   26.93    1.000
# # 1 DF_RSuccDate(x)          100 1196.95   44.447
#
# # benchmarking CSuccDate and RSuccDate with a data frame with 100 rows
# #              test replications elapsed relative
# # 2 DF_CSuccDate(y)          100    2.64    1.000
# # 1 DF_RSuccDate(x)          100  120.65   45.701
#
######   CSuccDate is about 45-times faster than RSuccDate   ###########################################################################
######                                                       ###########################################################################
########################################################################################################################################
#
#
#
#_______________________________________________________________________________________________________________________________________
#***************************************************************************************************************************************
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
########################################################################################################################################
######------>>>                                                               <<<------#################################################
####------>>>   comparing two methods that find the intersection of two sets    <<<------###############################################
######------>>>                                                               <<<------#################################################
#
# # method 1 is using intersect() of {base}
#
# USINGintersect<-function(x) {
#   for (i in c(1:nrow(x))) {
#     AnnivDates<-c(na.omit(as.numeric(x[i,c(9:29)])))
#     if (length(AnnivDates)>5) {
#       FIPD<-as.numeric(AnnivDates[2])
#       LIPD<-as.numeric(AnnivDates[5])
#     } else {
#       FIPD<-as.numeric(AnnivDates[1])
#       LIPD<-as.numeric(AnnivDates[1])
#     }
#     interTest<-intersect(which(!(AnnivDates<FIPD)),which(!(AnnivDates>LIPD)))
#     for (j in c(1:length(interTest))) {
#       x[i,j+29]<-interTest[j]
#     }
#   }
#   return(x)
# }
#
# # method 2 is using which() of {base}
#
# USINGwhich<-function(x) {
#   for (i in c(1:nrow(x))) {
#     AnnivDates<-c(na.omit(as.numeric(x[i,c(9:29)])))
#     if (length(AnnivDates)>5) {
#       FIPD<-as.numeric(AnnivDates[2])
#       LIPD<-as.numeric(AnnivDates[5])
#     } else {
#       FIPD<-as.numeric(AnnivDates[1])
#       LIPD<-as.numeric(AnnivDates[1])
#     }
#     whichTest<-which((!(AnnivDates<FIPD))&(!(AnnivDates>LIPD)))
#     for (j in c(1:length(whichTest))) {
#       x[i,j+29]<-whichTest[j]
#     }
#   }
#   return(x)
# }
#
# # creating a test data frame
# years1<-sample(c(1990:1999),10000,replace = TRUE)
# months1<-sample(c(1:12),10000,replace = TRUE)
# days1<-sample(c(1:28),10000,replace = TRUE)
# years2<-sample(c(2000:2010),10000,replace = TRUE)
# months2<-sample(c(1:12),10000,replace = TRUE)
# days2<-sample(c(1:28),10000,replace = TRUE)
# INTERSECTtest<-data.frame(years1,months1,days1,years2,months2,days2)
# INTERSECTtest$Date1<-as.Date(gsub("[ //s]","",paste(INTERSECTtest$years1,"-",INTERSECTtest$months1,"-",INTERSECTtest$days1)))
# INTERSECTtest$Date2<-as.Date(gsub("[ //s]","",paste(INTERSECTtest$years2,"-",INTERSECTtest$months2,"-",INTERSECTtest$days2)))
# for (i in c(1:nrow(INTERSECTtest))) {
#   helpVect<-c(seq(INTERSECTtest$Date1[i],INTERSECTtest$Date2[i],by="year"))
#   for (j in c(1:length(helpVect))) {
#     INTERSECTtest[i,j+8]<-as.Date(as.numeric(helpVect[j]),origin = "1970-1-1")
#   }
# }
# for (i in c(9:ncol(INTERSECTtest))) {
#   INTERSECTtest[,i]<-as.Date(INTERSECTtest[,i],origin = "1970-1-1")
# }
#
#
# # checking whether both functions do exactly the same
# DF1<-USINGintersect(INTERSECTtest)
# DF2<-USINGwhich(INTERSECTtest)
#
# length(which(DF1[,30]!=DF2[,30]))
# length(which(DF1[,31]!=DF2[,31]))
# length(which(DF1[,32]!=DF2[,32]))
# length(which(DF1[,33]!=DF2[,33]))
#
# library(rbenchmark)
#
# benchmark(USINGintersect(INTERSECTtest), USINGwhich(INTERSECTtest), order="relative")[,1:4]
#
# # benchmarking USINGintersect and USINGwhich with a data frame with 1000 rows
# #                            test replications elapsed relative
# # 2     USINGwhich(INTERSECTtest)          100   60.81    1.000
# # 1 USINGintersect(INTERSECTtest)          100   62.29    1.024
#
# # benchmarking timeDateAtomGet and StringAtomGet with a data frame with 10000 rows
# #                            test replications elapsed relative
# # 2     USINGwhich(INTERSECTtest)          100 1338.98    1.000
# # 1 USINGintersect(INTERSECTtest)          100 1340.44    1.001
#
######   USINGwhich is marginally faster than USINGintersect   #########################################################################
######        the advantage of USINGwhich decreases            #########################################################################
######              with increaseing data size                 #########################################################################
######                                                         #########################################################################
########################################################################################################################################
#
#
#
#_______________________________________________________________________________________________________________________________________
#***************************************************************************************************************************************
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
########################################################################################################################################
######------>>>                                                                             <<<------###################################
####------>>>   comparing two methods that extract the atoms of an object of class "Date"     <<<------#################################
######------>>>                                                                             <<<------###################################
#
# # method 1 is part of the package timeDate
#
# library(timeDate)
#
# timeDateAtomGet<-function(x) {
#   for (i in c(1:nrow(x))) {
#     x$Atom1[i]<-atoms(as.timeDate(x$Date[i]))[1]
#     x$Atom2[i]<-atoms(as.timeDate(x$Date[i]))[2]
#     x$Atom3[i]<-atoms(as.timeDate(x$Date[i]))[3]
#   }
#   return(x)
# }
#
# # method 2 uses string operations of base
# StringAtomGet<-function(x) {
#   for (i in c(1:nrow(x))) {
#     AtomVector<-as.numeric(unlist(strsplit(as.character(x$Date[i]),split = "-")))
#     x$Atom1[i]<-AtomVector[1]
#     x$Atom2[i]<-AtomVector[2]
#     x$Atom3[i]<-AtomVector[3]
#   }
#   return(x)
# }
#
# # creating a test data frame
# years<-sample(c(1000:3000),1000,replace = TRUE)
# months<-sample(c(1:12),1000,replace = TRUE)
# days<-sample(c(1:28),1000,replace = TRUE)
# AtomGetTest<-data.frame(years,months,days)
# AtomGetTest$Date<-as.Date(gsub("[ //s]","",paste(AtomGetTest$years,"-",AtomGetTest$months,"-",AtomGetTest$days)))
#
# # checking whether both functions do exactly the same
# DF1<-timeDateAtomGet(AtomGetTest)
# DF2<-StringAtomGet(AtomGetTest)
#
# length(which(DF1$Atom1!=DF2$Atom1))
# length(which(DF1$Atom2!=DF2$Atom2))
# length(which(DF1$Atom3!=DF2$Atom3))
#
# library(rbenchmark)
#
# # benchmarking the functions
# benchmark(timeDateAtomGet(AtomGetTest), StringAtomGet(AtomGetTest), order="relative")[,1:4]
#
# # benchmarking timeDateAtomGet and StringAtomGet with a data frame with 1000 rows
# #                           test replications elapsed relative
# # 2   StringAtomGet(AtomGetTest)          100   15.80    1.000
# # 1 timeDateAtomGet(AtomGetTest)          100  465.29   29.449
#
# # benchmarking timeDateAtomGet and StringAtomGet with a data frame with 5000 rows
# #                           test replications elapsed relative
# # 2   StringAtomGet(AtomGetTest)          100  108.75    1.000
# # 1 timeDateAtomGet(AtomGetTest)          100 2399.40   22.063
#
######   StringAtomGet is significantly faster than timeDateAtomGet   ##################################################################
######           the advantage of StringAtomGet decreases             ##################################################################
######                 with increaseing data size                     ##################################################################
######                                                                ##################################################################
########################################################################################################################################
#
#
#
#_______________________________________________________________________________________________________________________________________
#***************************************************************************************************************************************
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
########################################################################################################################################
#######------>>>                                       <<<------########################################################################
####------->>>   testing and benchmarking CppPrevDate    <<<------######################################################################
#######------>>>                                       <<<------########################################################################
#
#
#
#################################################################################################################################
######                                                                    #######################################################
######   testing whether RPrevDate and CppPrevDate do exactly the same    #######################################################
#
# library(timeDate)
#
# years1<-sample(c(1000:3000),1000,replace = TRUE)
# years2<-sample(c(1000:3000),1000,replace = TRUE)
# years3<-sample(c(1000:3000),1000,replace = TRUE)
# months1<-sample(c(1:12),1000,replace = TRUE)
# months2<-sample(c(1:12),1000,replace = TRUE)
# months3<-sample(c(1:12),1000,replace = TRUE)
# days1<-sample(c(1:28),1000,replace = TRUE)
# days2<-sample(c(1:28),1000,replace = TRUE)
# days3<-sample(c(1:28),1000,replace = TRUE)
# CpY<-sample(c(1,2,3,4,6,12),1000,replace = TRUE)
# PrevDateTest<-data.frame(years1,months1,days1,years2,months2,days2,years3,months3,days3,CpY)
# PrevDateTest$Date1<-as.Date(gsub("[ //s]","",paste(PrevDateTest$years1,"-",PrevDateTest$months1,"-",PrevDateTest$days1)))
# PrevDateTest$Date2<-as.Date(gsub("[ //s]","",paste(PrevDateTest$years2,"-",PrevDateTest$months2,"-",PrevDateTest$days2)))
# PrevDateTest$Date3<-as.Date(gsub("[ //s]","",paste(PrevDateTest$years3,"-",PrevDateTest$months3,"-",PrevDateTest$days3)))
#
# # creating the anniversary date preceding Em in R
# RPrevDate<-function(AnnivDates1,Em,Mat,CpY) {
#   gap<-round(abs(as.numeric(difftime(as.Date(ISOdatetime(t(atoms(as.timeDate(Em))[1]),t(atoms(as.timeDate(Em))[2]),15,12,0,0)),
#                                      as.Date(ISOdatetime(t(atoms(as.timeDate(Em))[1]+1),t(atoms(as.timeDate(Em))[2]),15,12,0,0)))))/CpY,digits=0)
#   if (!(AnnivDates1<Em)) {
#     PrevDate<-as.Date(ISOdatetime(t(atoms(as.timeDate(AnnivDates1))[1]),t(atoms(as.timeDate(AnnivDates1))[2]),15,12,0,0))
#     PrevDate<-PrevDate-gap
#     if ((as.Date(timeLastDayInMonth(Mat))==Mat)==TRUE) {
#       PrevDate<-as.Date(timeLastDayInMonth(PrevDate))
#     } else {
#       PrevDate_A<-as.Date(ISOdatetime(t(atoms(as.timeDate(PrevDate))[1]),t(atoms(as.timeDate(PrevDate))[2]),t(atoms(as.timeDate(Mat))[3]),12,0,0))
#       if (is.na(PrevDate_A)) {
#         PrevDate<-as.Date(timeLastDayInMonth(PrevDate))
#       } else {
#         PrevDate<-PrevDate_A
#       }
#     }
#   } else {
#     PrevDate<-NA
#   }
#   PrevDate<-as.Date(PrevDate)
#   return(PrevDate)
# }
#
# # creating the anniversary date preceding Em in C++
# CPrevDate<-function(x) {
#   out<-CppPrevDate(x)
#   if (any(is.na(out))) {
#     out<-as.Date(NA)
#   } else {
#     out<-as.Date(paste(out,collapse="-"),origin="1970-01-01")
#   }
#   return(out)
# }
#
# for (i in c(1:nrow(PrevDateTest))) {
#   PrevDateTest$RPrevDate[i]<-RPrevDate(PrevDateTest$Date1[i],PrevDateTest$Date2[i],PrevDateTest$Date3[i],PrevDateTest$CpY[i])
# }
#
# for (i in c(1:nrow(PrevDateTest))) {
#   PrevDateTest$CPrevDate[i]<-CPrevDate(c(PrevDateTest$years1[i],PrevDateTest$months1[i],PrevDateTest$days1[i],
#                                 PrevDateTest$years2[i],PrevDateTest$months2[i],PrevDateTest$days2[i],
#                                 PrevDateTest$years3[i],PrevDateTest$months3[i],PrevDateTest$days3[i],
#                                 PrevDateTest$CpY[i]))
# }
#
# for (i in c(1:nrow(PrevDateTest))) {
#   if (all(is.na(PrevDateTest$CPrevDate[i]),is.na(PrevDateTest$RPrevDate[i]))) {
#     PrevDateTest$Check[i]<-1
#   } else {
#     if (PrevDateTest$CPrevDate[i]==PrevDateTest$RPrevDate[i]) {
#       PrevDateTest$Check[i]<-1
#     } else {
#       PrevDateTest$Check[i]<-0
#     }
#   }
# }
#
# length(which(PrevDateTest$Check==0))
#
#
######   HURRRAAAAYYYY!!! RPrevDate and CPrevDate do exactly the same!    #######################################################
######                                                                    #######################################################
#################################################################################################################################
###^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^###
###vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv###
#################################################################################################################################
######                                            ###############################################################################
######   benchmarking RPrevDate and CPrevDate     ###############################################################################
#
# years1<-sample(c(1000:3000),100,replace = TRUE)
# years2<-sample(c(1000:3000),100,replace = TRUE)
# years3<-sample(c(1000:3000),100,replace = TRUE)
# months1<-sample(c(1:12),100,replace = TRUE)
# months2<-sample(c(1:12),100,replace = TRUE)
# months3<-sample(c(1:12),100,replace = TRUE)
# days1<-sample(c(1:28),100,replace = TRUE)
# days2<-sample(c(1:28),100,replace = TRUE)
# days3<-sample(c(1:28),100,replace = TRUE)
# CpY<-sample(c(1,2,3,4,6,12),100,replace = TRUE)
# PrevDateTest<-data.frame(years1,months1,days1,years2,months2,days2,years3,months3,days3,CpY)
# PrevDateTest$Date1<-as.Date(gsub("[ //s]","",paste(PrevDateTest$years1,"-",PrevDateTest$months1,"-",PrevDateTest$days1)))
# PrevDateTest$Date2<-as.Date(gsub("[ //s]","",paste(PrevDateTest$years2,"-",PrevDateTest$months2,"-",PrevDateTest$days2)))
# PrevDateTest$Date3<-as.Date(gsub("[ //s]","",paste(PrevDateTest$years3,"-",PrevDateTest$months3,"-",PrevDateTest$days3)))
#
#
# library(rbenchmark)
#
# DF_RPrevDate<-function(x) {
#   OUT<-RPrevDate(x[1,1],x[1,2],x[1,3],x[1,4])
#   for (i in c(2:nrow(x))) {
#     OUTnext<-RPrevDate(x[i,1],x[i,2],x[i,3],x[i,4])
#     OUT<-append(OUT,OUTnext)
#   }
#   return(OUT)
# }
#
# DF_CPrevDate<-function(x) {
#   OUT<-CPrevDate(c(x[1,1],x[1,2],x[1,3],x[1,4],x[1,5],x[1,6],x[1,7],x[1,8],x[1,9],x[1,10]))
#   for (i in c(2:nrow(x))) {
#     OUTnext<-CPrevDate(c(x[i,1],x[i,2],x[i,3],x[i,4],x[i,5],x[i,6],x[i,7],x[i,8],x[i,9],x[i,10]))
#     OUT<-append(OUT,OUTnext)
#   }
#   return(OUT)
# }
#
# x<-PrevDateTest[,c(11:13,10)]
# y<-PrevDateTest[,c(1:10)]
#
# # check that DF_RPrevDate is the same as DF_CPrevDate
# stopifnot(all.equal(DF_RPrevDate(x), DF_CPrevDate(y)))
#
# benchmark(DF_RPrevDate(x), DF_CPrevDate(y), order="relative")[,1:4]
#
# # benchmarking CPrevDate and RPrevDate with a data frame with 1000 rows
# #              test replications elapsed relative
# # 2 DF_CPrevDate(y)          100   26.94    1.000
# # 1 DF_RPrevDate(x)          100 1201.85   44.612
#
# # benchmarking CPrevDate and RPrevDate with a data frame with 100 rows
# #              test replications elapsed relative
# # 2 DF_CPrevDate(y)          100    2.56    1.000
# # 1 DF_RPrevDate(x)          100  113.39   44.293
#
######   CPrevDate is about 45-times faster than RPrevDate   ###########################################################################
######                                                       ###########################################################################
########################################################################################################################################
#
#
#
#_______________________________________________________________________________________________________________________________________
#***************************************************************************************************************************************
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
########################################################################################################################################
######------>>>                                    <<<------############################################################################
####------>>>   testing and benchmarking NumToDate   <<<------##########################################################################
######------>>>                                    <<<------############################################################################
#
#
#
#################################################################################################################################
######                                                                ###########################################################
######   testing whether NumToDate and as.Date do exactly the same    ###########################################################
#
# x<-runif(1000)
# y<-runif(1000)
# x<-trunc(x*1000000)
# y<-trunc(y*1000000)
# y<-(-1)*y
# x<-c(x,y)
# originYear<-1099
# originMonth<-7
# originDay<-15
#
# MyNumToDate<-function(x) {
#   DateVector<-c(NumToDate(c(x[1],originYear,originMonth,originDay)))
#   DateVector<-as.character(gsub("[ \\s]","",paste(DateVector[1],"-",DateVector[2],"-",DateVector[3])))
#   for (i in c(2:length(x))) {
#     DateVectorNext<-c(NumToDate(c(x[i],originYear,originMonth,originDay)))
#     DateVectorNext<-as.character(gsub("[ \\s]","",paste(DateVectorNext[1],"-",DateVectorNext[2],"-",DateVectorNext[3])))
#     DateVector<-append(DateVector,DateVectorNext)
#   }
#   return(DateVector)
# }
#
#
# RNumToDate<-function(x) {
#   DateVector<-as.Date(x[1],origin=gsub("[ //s]","",paste(originYear,"-",originMonth,"-",originDay)))
#   DateVector<-unlist(strsplit(as.character(DateVector),""))
#   if (DateVector[1]=="-") {
#     NegYear<-1
#     DateVector<-DateVector[-1]
#   } else {
#     NegYear<-0
#   }
#   DateVector<-paste(DateVector, collapse="")
#   DateVector<-as.numeric(unlist(strsplit(as.character(DateVector),"-")))
#   if (NegYear==1) {
#     DateVector<-as.character(gsub("[ \\s]","",paste("-",DateVector[1],"-",DateVector[2],"-",DateVector[3])))
#   } else {
#     DateVector<-as.character(gsub("[ \\s]","",paste(DateVector[1],"-",DateVector[2],"-",DateVector[3])))
#   }
#   for (i in c(2:length(x))) {
#     DateVectorNext<-as.Date(x[i],origin=gsub("[ //s]","",paste(originYear,"-",originMonth,"-",originDay)))
#     DateVectorNext<-unlist(strsplit(as.character(DateVectorNext),""))
#     if (DateVectorNext[1]=="-") {
#       NegYear<-1
#       DateVectorNext<-DateVectorNext[-1]
#     } else {
#       NegYear<-0
#     }
#     DateVectorNext<-paste(DateVectorNext, collapse="")
#     DateVectorNext<-as.numeric(unlist(strsplit(as.character(DateVectorNext),"-")))
#     if (NegYear==1) {
#       DateVectorNext<-as.character(gsub("[ \\s]","",paste("-",DateVectorNext[1],"-",DateVectorNext[2],"-",DateVectorNext[3])))
#     } else {
#       DateVectorNext<-as.character(gsub("[ \\s]","",paste(DateVectorNext[1],"-",DateVectorNext[2],"-",DateVectorNext[3])))
#     }
#     DateVector<-append(DateVector,DateVectorNext)
#   }
#   return(DateVector)
# }
#
# MY<-MyNumToDate(x)
# RsFunc<-RNumToDate(x)
# Check<-data.frame(MY,RsFunc)
# Check[,1]<-as.character(Check[,1])
# Check[,2]<-as.character(Check[,2])
#
# for (i in c(1:nrow(Check))) {
#   if (is.na(Check[i,1])) {
#     Check[i,3]<-0
#   } else {
#     if (Check[i,1]==Check[i,2]) {
#       Check[i,3]<-1
#     } else {
#       Check[i,3]<-0
#     }
#   }
# }
#
# length(which(Check[,3]==0))
#
# # check that MyNumToDate is the same as RNumToDate
# stopifnot(all.equal(MyNumToDate(x), RNumToDate(x)))
#
#
######   HURRRAAAAYYYY!!! NumToDate and as.Date do exactly the same!    #########################################################
######                                                                  #########################################################
#################################################################################################################################
###^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^###
###vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv###
#################################################################################################################################
######                                          #################################################################################
######   benchmarking NumToDate and as.Date     #################################################################################
#
# MyNumToDate<-function(x) {
#   DateVector<-c(NumToDate(c(x[1],originYear,originMonth,originDay)))
#   for (i in c(2:length(x))) {
#     DateVectorNext<-c(NumToDate(c(x[i],originYear,originMonth,originDay)))
#     DateVector<-append(DateVector,DateVectorNext)
#   }
#   return(DateVector)
# }
#
# RNumToDate<-function(x) {
#   DateVector<-c(as.Date(x[1],origin=gsub("[ //s]","",paste(originYear,"-",originMonth,"-",originDay))))
#   for (i in c(2:length(x))) {
#     DateVectorNext<-as.Date(x[i],origin=gsub("[ //s]","",paste(originYear,"-",originMonth,"-",originDay)))
#     DateVector<-append(DateVector,DateVectorNext)
#   }
#   return(DateVector)
# }
#
#
# library(rbenchmark)
# #
# #
# # # benchmarking MyNumToDate and RNumToDate with dates AD and dates BC
# x<-runif(1000)
# y<-runif(1000)
# x<-trunc(x*1000000)
# y<-trunc(y*1000000)
# y<-(-1)*y
# x<-c(x,y)
# originYear<-1985
# originMonth<-4
# originDay<-21
# benchmark(MyNumToDate(x), RNumToDate(x), order="relative")[,1:4]
#
# #             test replications elapsed relative
# # 1 MyNumToDate(x)          100    5.71     1.00
# # 2  RNumToDate(x)          100   34.89     6.11
#
#
#
######   NumToDate is about 6-times faster than as.Date   ##############################################################################
######                                                    ##############################################################################
########################################################################################################################################
#
#
#
#_______________________________________________________________________________________________________________________________________
#***************************************************************************************************************************************
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
########################################################################################################################################
######------>>>                 <<<------###############################################################################################
####------>>>   testing DayDiff   <<<------#############################################################################################
######------>>>                 <<<------###############################################################################################
#
#
#
#################################################################################################################################
######                                                                ###########################################################
######   testing whether DayDiff and difftime do exactly the same     ###########################################################
#
# years1<-trunc(runif(1000)*10000)
# years2<-trunc(runif(1000)*10000)
# months1<-sample(c(1:12),1000,replace = TRUE)
# months2<-sample(c(1:12),1000,replace = TRUE)
# days1<-sample(c(1:28),1000,replace = TRUE)
# days2<-sample(c(1:28),1000,replace = TRUE)
# DayDiffTest<-data.frame(years1,months1,days1,years2,months2,days2)
# DayDiffTest$Date1<-as.Date(gsub("[ //s]","",paste(DayDiffTest$years1,"-",DayDiffTest$months1,"-",DayDiffTest$days1)))
# DayDiffTest$Date2<-as.Date(gsub("[ //s]","",paste(DayDiffTest$years2,"-",DayDiffTest$months2,"-",DayDiffTest$days2)))
# DayDiffTest$RFunc<-as.numeric(difftime(DayDiffTest$Date2,DayDiffTest$Date1))
#
# for (i in c(1:nrow(DayDiffTest))) {
#   DayDiffTest$MyFunc[i]<-DayDiff(c(DayDiffTest$years1[i],DayDiffTest$months1[i],DayDiffTest$days1[i],
#                                 DayDiffTest$years2[i],DayDiffTest$months2[i],DayDiffTest$days2[i]))
# }
#
# for (i in c(1:nrow(DayDiffTest))) {
#   if (is.na(DayDiffTest$MyFunc[i])) {
#     DayDiffTest$Check[i]<-0
#   } else {
#     if (DayDiffTest$RFunc[i]==DayDiffTest$MyFunc[i]) {
#       DayDiffTest$Check[i]<-1
#     } else {
#       DayDiffTest$Check[i]<-0
#     }
#   }
# }
#
# length(which(DayDiffTest$Check==0))
#
######    HURRRAAAAYYYY!!! DayDiff and difftime do exactly the same!    #########################################################
######                                                                  #########################################################
#################################################################################################################################
###^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^###
###vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv###
#################################################################################################################################
######                                          #################################################################################
######   benchmarking DayDiff and difftime      #################################################################################
#
# library(rbenchmark)
#
# MyDayDiff<-function(x) {
#   MyDayDiffTest<-DayDiff(c(x[1,1],x[1,2],x[1,3],x[1,4],x[1,5],x[1,6]))
#   for (i in c(2:nrow(x))) {
#     MyDayDiffTestNext<-DayDiff(c(x[i,1],x[i,2],x[i,3],x[i,4],x[i,5],x[i,6]))
#     MyDayDiffTest<-append(MyDayDiffTest,MyDayDiffTestNext)
#   }
#   return(MyDayDiffTest)
# }
#
# Rdifftime<-function(x) {
#   RdifftimeTest<-as.numeric(difftime(x[1,2],x[1,1]))
#   for (i in c(2:nrow(x))) {
#     RdifftimeTestNext<-as.numeric(difftime(x[i,2],x[i,1]))
#     RdifftimeTest<-append(RdifftimeTest,RdifftimeTestNext)
#   }
#   return(RdifftimeTest)
# }
#
# x<-DayDiffTest[,c(1:6)]
# y<-DayDiffTest[,c(7,8)]
#
# # check that MyNumToDate is the same as RNumToDate
# stopifnot(all.equal(MyDayDiff(x), Rdifftime(y)))
#
#
# benchmark(MyDayDiff(x), Rdifftime(y), order="relative")[,1:4]
#
# #           test replications elapsed relative
# # 1 MyDayDiff(x)          100   11.39    1.000
# # 2 Rdifftime(y)          100   13.98    1.227
#
######   DayDiff is about 1.2-times faster than difftime  ##############################################################################
######                                                    ##############################################################################
########################################################################################################################################
#
#
#
#_______________________________________________________________________________________________________________________________________
#***************************************************************************************************************************************
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
########################################################################################################################################
######------>>>                                                  <<<------##############################################################
####------>>>  testing outdated DayCountCon and outdated PayCalc   <<<------############################################################
######------>>>                                                  <<<------##############################################################
#
# AtomVec<-c(2010,8,31,2010,11,2,2011,2,28,2,1,2015,8,31,100,0.15)
#
# DayCountCon(AtomVec)
#
# DayDiff(c(2010,8,31,2011,1,1))
#
# CpY<-4
# DCC<-8
# AtomVec1<-c(2010,2,28,2010,4,21,2010,5,31,CpY,DCC,2015,8,31,100,0.15)
# AtomVec2<-c(2010,5,31,2010,7,27,2010,8,31,CpY,DCC,2015,8,31,100,0.15)
# AtomVec3<-c(2010,8,31,2011,11,2,2010,11,30,CpY,DCC,2015,8,31,100,0.15)
# AtomVec4<-c(2010,11,30,2011,12,24,2011,2,28,CpY,DCC,2015,8,31,100,0.15)
#
# DF<-data.frame(t(AtomVec1))
# DF<-rbind(DF,AtomVec2,AtomVec3,AtomVec4)
# DF<-as.matrix(DF)
#
# PayCalc(DF)
#
#_______________________________________________________________________________________________________________________________________
#***************************************************************************************************************************************
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
########################################################################################################################################
######------>>>                                                  <<<------##############################################################
####------>>>   writing a method to use PayCalc in an R Function   <<<------############################################################
######------>>>                                                  <<<------##############################################################
#
# Em<-as.Date("2009-6-15")
# CD1<-as.Date("2010-2-28")
# CD2<-as.Date("2010-5-31")
# CD3<-as.Date("2010-8-31")
# CD4<-as.Date("2010-11-30")
# Mat<-as.Date("2010-12-24")
# CpY<-4
# Coup<-0.15
# RV<-100
# DCC<-1
#
# RealDates<-c(Em,CD1,CD2,CD3,CD4,Mat)
#
# D1<-RealDates[-length(RealDates)]
# D2<-RealDates[-1]
# D3<-RealDates[-1]
# MatVec<-rep(Mat,length(RealDates)-1)
# CpYVec<-rep(CpY,length(RealDates)-1)
# CoupVec<-rep(Coup,length(RealDates)-1)
# RVVec<-rep(RV,length(RealDates)-1)
# DCCVec<-rep(DCC,length(RealDates)-1)
# D1Atoms<-as.numeric(unlist(strsplit(as.character(D1),split="-")))
# D2Atoms<-as.numeric(unlist(strsplit(as.character(D2),split="-")))
# D3Atoms<-as.numeric(unlist(strsplit(as.character(D3),split="-")))
# MatAtoms<-as.numeric(unlist(strsplit(as.character(MatVec),split="-")))
# D1Matrix<-matrix(t(D1Atoms),nrow=length(RealDates)-1,ncol=3,byrow=TRUE)
# D2Matrix<-matrix(t(D2Atoms),nrow=length(RealDates)-1,ncol=3,byrow=TRUE)
# D3Matrix<-matrix(t(D3Atoms),nrow=length(RealDates)-1,ncol=3,byrow=TRUE)
# MatMatrix<-matrix(t(MatAtoms),nrow=length(RealDates)-1,ncol=3,byrow=TRUE)
# CoupSchedMatrix<-cbind(D1Matrix,D2Matrix,D3Matrix,CpY,DCC,MatMatrix,RV,Coup)
# CoupPayments<-PayCalc(CoupSchedMatrix)
#
# # this is the input required by PayCalc: first row: Y1,M1,D1,Y2,M2,D2,Y3,M3,D3,CpY,DCC,YMat,MMat,DMat,RV,Coup
#
##############   the transformations above have to be done within the     ##############################################################
################               the functon AnnivDates                 ##################################################################
######################                                          ########################################################################
########################################################################################################################################



*/











