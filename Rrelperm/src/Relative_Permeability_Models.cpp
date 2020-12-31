
#include "RcppArmadillo.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]


//' @importFrom Rcpp evalCpp
//' @importFrom Rdpack reprompt
//' @useDynLib Rrelperm, .registration=TRUE


//******************************************************************************
// 2-PHASE BROOKS-COREY, KROW


// [[Rcpp::export]]
NumericMatrix kr2p_ow_cpp(double SWCON, double SWCRIT, double SOIRW, double SORW,
                             double KRWIRO, double KROCW, double NW, double NOW, int NP) {

   double threshold = 1e-6;
   int len = NP;
   if (len > 501) {
      len = 501;
   }
   arma::mat mat(len,4);
   arma::vec SW = arma::regspace(0., 1. / (NP - 1), 1.);
   arma::vec SO(len);
   for (int i = 0; i < len; i++) {
      mat(i,0) = SW(i);
      SO(i) = 1 - SW(i);
      mat(i,1) = SO(i);
      if (SW(i) <= SWCRIT) {
         mat(i,2) = 0;
      } else if (SW(i) >= (1 - SOIRW)) {
         mat(i,2) = KRWIRO;
      } else {
         mat(i,2) = KRWIRO * std::pow(((SW(i) - SWCRIT) / (1 - SWCRIT - SOIRW)), NW);
         if (std::abs(mat(i,2)) < threshold) {
            mat(i,2) = 0.0;
         }
      }
      if (SO(i) <= SORW) {
         mat(i,3) = 0;
      } else if (SO(i) >= (1 - SWCON)) {
         mat(i,3) = KROCW;
      } else {
         mat(i,3) = KROCW * std::pow(((SO(i) - SORW) / (1 - SWCON - SORW)), NOW);
         if (std::abs(mat(i,3)) < threshold) {
            mat(i,3) = 0.0;
         }
      }
   }
   NumericMatrix mat22 = wrap(mat);
   CharacterVector colname(4);
   colname = {"Sw", "So", "Krw", "Krow"};
   colnames(mat22) = colname;
   return(mat22);
}


//******************************************************************************
// 2-PHASE BROOKS-COREY, KRGL


// [[Rcpp::export]]
NumericMatrix kr2p_gl_cpp(double SWCON, double SOIRG, double SORG, double SGCON, double SGCRIT,
                  double KRGCL, double KROGCG, double NG, double NOG, int NP) {

   double threshold = 1e-6;
   int len = NP;
   if (len > 501) {
      len = 501;
   }
   arma::mat mat(len,4);
   arma::vec SG = arma::regspace(0., 1.0 / (NP - 1), 1.);
   arma::vec SL(len);
   double SLCON = SWCON + SOIRG;
   double SLRG = SWCON + SORG;
   for (int i = 0; i < len; i++) {
      mat(i,0) = SG(i);
      SL(i) = 1 - SG(i);
      mat(i,1) = SL(i);
      if (SG(i) <= SGCRIT) {
         mat(i,2) = 0;
      } else if (SG(i) >= (1 - SLCON)) {
         mat(i,2) = KRGCL;
      } else {
         mat(i,2) = KRGCL * std::pow(((SG(i) - SGCRIT) / (1 - SGCRIT - SLCON)), NG);
         if (std::abs(mat(i,2)) < threshold) {
            mat(i,2) = 0.0;
         }
      }
      if (SL(i) <= SLRG) {
         mat(i,3) = 0;
      } else if (SL(i) >= (1 - SGCON)) {
         mat(i,3) = KROGCG;
      } else {
         mat(i,3) = KROGCG * std::pow(((SL(i) - SLRG) / (1 - SGCON - SLRG)), NOG);
         if (std::abs(mat(i,3)) < threshold) {
            mat(i,3) = 0.0;
         }
      }
   }
   NumericMatrix mat22 = wrap(mat);
   CharacterVector colname(4);
   colname = {"Sg", "Sl", "Krg", "Krog"};
   colnames(mat22) = colname;
   return(mat22);
}




//******************************************************************************
// 3-PHASE STONE I, SO


// [[Rcpp::export]]
NumericMatrix kr3p_StoneI_So_cpp(double SWCON, double SWCRIT, double SOIRW, double SORW,
                          double SOIRG, double SORG, double SGCON, double SGCRIT,
                          double KRWIRO, double KROCW, double KRGCL,
                          double NW, double NOW, double NG, double NOG, int NP) {

   const double threshold = 1e-6;
   int len = NP;
   if (len > 501) {
      len = 501;
   }
   double KROGCG = KROCW;
   arma::vec SW = arma::regspace(0., 1. / (len - 1), 1.);
   arma::vec SG = arma::regspace(0., 1. / (len - 1), 1.);
   double SLRG = SWCON + SORG;
   int lent = len * (len + 1) / 2;
   int counter_1 = 0;
   int counter_2 = 0;
   arma::vec SW3(lent);
   arma::vec SG3(lent);
   arma::vec SO3(lent);
   arma::mat mat2(lent,4);
   double temp, KROW, KROG;
   double a, SOm, SOstar, SWstar, SGstar;
   for (int i = 0; i < len; i++) {
      for (int j = 0; j < (len - i); j++) {
         SW3(counter_2 + j) = SW(i);
         SG3(counter_2 + j) = SG(j);
         SO3(counter_2 + j) = 1 - SW(i) - SG(j);
         mat2(counter_2 + j,0) = SW3(counter_2 + j);
         mat2(counter_2 + j,1) = SG3(counter_2 + j);
         mat2(counter_2 + j,2) = SO3(counter_2 + j);
         a = SG3(counter_2 + j) / (1 - SWCON - SORG);
         SOm = (1. - a) * SORW + a * SORG;
         SOstar = (SO3(counter_2 + j) - SOm) / (1 - SWCON - SOm);
         SWstar = (SW3(counter_2 + j) - SWCON) / (1 - SWCON - SOm);
         SGstar = SG3(counter_2 + j) / (1 - SWCON - SOm);
         temp = SO3(counter_2 + j);
         if (temp <= SORW) {
            KROW = 0.;
         } else if (temp >= (1. - SWCON)) {
            KROW = KROCW;
         } else {
            KROW = KROCW * std::pow(((temp - SORW) / (1. - SWCON - SORW)), NOW);
            if (std::abs(KROW) < threshold) {
               KROW = 0.0;
            }
         }
         temp = SO3(counter_2 + j) + SWCON;
         if (temp <= SLRG) {
            KROG = 0.;
         } else if (temp >= (1. - SGCON)) {
            KROG = KROGCG;
         } else {
            KROG = KROGCG * std::pow(((temp - SLRG) / (1. - SGCON - SLRG)), NOG);
            if (std::abs(KROG) < threshold) {
               KROG = 0.0;
            }
         }
         temp = SOstar * KROW * KROG;
         if (temp <= 0.) {
            mat2(counter_2 + j,3) = 0.;
         } else {
            mat2(counter_2 + j,3) = SOstar * KROW * KROG / (KROCW * (1. - SWstar) * (1. - SGstar));
         }
         if (mat2(counter_2 + j,3) > 1.) {
            mat2(counter_2 + j,3) = 1.0;
         }
         if (mat2(counter_2 + j,3) < 0.) {
            mat2(counter_2 + j,3) = 0.0;
         }
         if (mat2(counter_2 + j,3) < threshold) {
            mat2(counter_2 + j,3) = 0.0;
         }
      }
      counter_1 = len - i;
      counter_2 = counter_2 + counter_1;
   }
   NumericMatrix mat22(lent,4);
   mat22 = Rcpp::wrap(mat2);
   CharacterVector colname(4);
   colname = {"Sw", "Sg", "So", "Kro"};
   colnames(mat22) = colname;
   return(mat22);
}


//******************************************************************************
// 3-PHASE STONE I, SWSG


// [[Rcpp::export]]
NumericMatrix kr3p_StoneI_SwSg_cpp(double SWCON, double SWCRIT, double SOIRW, double SORW,
                           double SOIRG, double SORG, double SGCON, double SGCRIT,
                           double KRWIRO, double KROCW, double KRGCL,
                           double NW, double NOW, double NG, double NOG, int NP) {

   const double threshold = 1e-6;
   int len = NP;
   if (len > 501) {
      len = 501;
   }
   double KROGCG = KROCW;
   arma::vec SW = arma::regspace(0., 1. / (len - 1), 1.);
   arma::vec SG = arma::regspace(0., 1. / (len - 1), 1.);
   double SLRG = SWCON + SORG;
   int lent = len * (len + 1) / 2;
   int counter_1 = 0;
   int counter_2 = 0;
   arma::vec SW3(lent);
   arma::vec SG3(lent);
   arma::vec SO3(lent);
   arma::mat mat2(lent,4);
   double temp, KROW, KROG;
   double a, SOm, SOstar, SWstar, SGstar;
   for (int i = 0; i < len; i++) {
      for (int j = 0; j < (len - i); j++) {
         SW3(counter_2 + j) = SW(i);
         SG3(counter_2 + j) = SG(j);
         SO3(counter_2 + j) = 1 - SW(i) - SG(j);
         mat2(counter_2 + j,0) = SW3(counter_2 + j);
         mat2(counter_2 + j,1) = SG3(counter_2 + j);
         mat2(counter_2 + j,2) = SO3(counter_2 + j);
         a = SG3(counter_2 + j) / (1 - SWCON - SORG);
         SOm = (1. - a) * SORW + a * SORG;
         SOstar = (SO3(counter_2 + j) - SOm) / (1 - SWCON - SOm);
         SWstar = (SW3(counter_2 + j) - SWCON) / (1 - SWCON - SOm);
         SGstar = SG3(counter_2 + j) / (1 - SWCON - SOm);
         temp = 1. - SW3(counter_2 + j) ;
         if (temp <= SORW) {
            KROW = 0.;
         } else if (temp >= (1 - SWCON)) {
            KROW = KROCW;
         } else {
            KROW = KROCW * std::pow(((temp - SORW) / (1 - SWCON - SORW)), NOW);
            if (std::abs(KROW) < threshold) {
               KROW = 0.0;
            }
         }
         temp = 1. - SG3(counter_2 + j);
         if (temp <= SLRG) {
            KROG = 0.;
         } else if (temp >= (1 - SGCON)) {
            KROG = KROGCG;
         } else {
            KROG = KROGCG * std::pow(((temp - SLRG) / (1 - SGCON - SLRG)), NOG);
            if (std::abs(KROG) < threshold) {
               KROG = 0.0;
            }
         }
         temp = SOstar * KROW * KROG;
         if (temp <= 0.) {
            mat2(counter_2 + j,3) = 0.;
         } else {
            mat2(counter_2 + j,3) = SOstar * KROW * KROG / (KROCW * (1 - SWstar) * (1 - SGstar));
         }
         if (mat2(counter_2 + j,3) > 1) {
            mat2(counter_2 + j,3) = 1.0;
         }
         if (mat2(counter_2 + j,3) < 0.) {
            mat2(counter_2 + j,3) = 0.0;
         }
         if (mat2(counter_2 + j,3) < threshold) {
            mat2(counter_2 + j,3) = 0.0;
         }
      }
      counter_1 = len - i;
      counter_2 = counter_2 + counter_1;
   }
   NumericMatrix mat22(lent,4);
   mat22 = Rcpp::wrap(mat2);
   CharacterVector colname(4);
   colname = {"Sw", "Sg", "So", "Kro"};
   colnames(mat22) = colname;
   return(mat22);
}



//******************************************************************************
// 3-PHASE STONE II, SO


// [[Rcpp::export]]
NumericMatrix kr3p_StoneII_So_cpp(double SWCON, double SWCRIT, double SOIRW, double SORW,
                          double SOIRG, double SORG, double SGCON, double SGCRIT,
                          double KRWIRO, double KROCW, double KRGCL,
                          double NW, double NOW, double NG, double NOG, int NP) {

   const double threshold = 1e-6;
   int len = NP;
   if (len > 501) {
      len = 501;
   }
   double KROGCG = KROCW;
   arma::vec SW = arma::regspace(0., 1. / (len - 1), 1.);
   arma::vec SG = arma::regspace(0., 1. / (len - 1), 1.);
   double SLCON = SWCON + SOIRG;
   double SLRG = SWCON + SORG;
   int lent = len * (len + 1) / 2;
   int counter_1 = 0;
   int counter_2 = 0;
   arma::vec SW3(lent);
   arma::vec SG3(lent);
   arma::vec SO3(lent);
   arma::mat mat2(lent,4);
   double temp, KRW, KRG, KROW, KROG;
   for (int i = 0; i < len; i++) {
      for (int j = 0; j < (len - i); j++) {
         SW3(counter_2 + j) = SW(i);
         SG3(counter_2 + j) = SG(j);
         SO3(counter_2 + j) = 1 - SW(i) - SG(j);
         mat2(counter_2 + j,0) = SW3(counter_2 + j);
         mat2(counter_2 + j,1) = SG3(counter_2 + j);
         mat2(counter_2 + j,2) = SO3(counter_2 + j);
         if (SW3(counter_2 + j)  <= SWCRIT) {
            KRW = 0.;
         } else if (SW3(counter_2 + j)  >= (1. - SOIRW)) {
            KRW = KRWIRO;
         } else {
            KRW = KRWIRO * std::pow(((SW3(counter_2 + j)  - SWCRIT) / (1. - SWCRIT - SOIRW)), NW);
            if (std::abs(KRW) < threshold) {
               KRW = 0.0;
            }
         }
         if (SG3(counter_2 + j) <= SGCRIT) {
            KRG = 0.;
         } else if (SG3(counter_2 + j) >= (1. - SLCON)) {
            KRG = KRGCL;
         } else {
            KRG = KRGCL * std::pow(((SG3(counter_2 + j) - SGCRIT) / (1. - SGCRIT - SLCON)), NG);
            if (std::abs(KRG) < threshold) {
               KRG = 0.0;
            }
         }
         temp = SO3(counter_2 + j);
         if (temp <= SORW) {
            KROW = 0.;
         } else if (temp >= (1. - SWCON)) {
            KROW = KROCW;
         } else {
            KROW = KROCW * std::pow(((temp - SORW) / (1. - SWCON - SORW)), NOW);
            if (std::abs(KROW) < threshold) {
               KROW = 0.0;
            }
         }
         temp = SO3(counter_2 + j) + SWCON;
         if (temp <= SLRG) {
            KROG = 0.;
         } else if (temp >= (1. - SGCON)) {
            KROG = KROGCG;
         } else {
            KROG = KROGCG * std::pow(((temp - SLRG) / (1. - SGCON - SLRG)), NOG);
            if (std::abs(KROG) < threshold) {
               KROG = 0.0;
            }
         }
         mat2(counter_2 + j,3) = KROCW * ((KROW / KROCW + KRW) * (KROG / KROCW + KRG) - KRW - KRG);
         if (mat2(counter_2 + j,3) < 0.) {
            mat2(counter_2 + j,3) = 0.;
         }
      }
      counter_1 = len - i;
      counter_2 = counter_2 + counter_1;
   }
   NumericMatrix mat22(lent,4);
   mat22 = Rcpp::wrap(mat2);
   CharacterVector colname(4);
   colname = {"Sw", "Sg", "So", "Kro"};
   colnames(mat22) = colname;
   return(mat22);
}



//******************************************************************************
// 3-PHASE STONE II, SWSG


// [[Rcpp::export]]
NumericMatrix kr3p_StoneII_SwSg_cpp(double SWCON, double SWCRIT, double SOIRW, double SORW,
                            double SOIRG, double SORG, double SGCON, double SGCRIT,
                            double KRWIRO, double KROCW, double KRGCL,
                            double NW, double NOW, double NG, double NOG, int NP) {

   const double threshold = 1e-6;
   int len = NP;
   if (len > 501) {
      len = 501;
   }
   double KROGCG = KROCW;
   arma::vec SW = arma::regspace(0., 1. / (len - 1), 1.);
   arma::vec SG = arma::regspace(0., 1. / (len - 1), 1.);
   double SLCON = SWCON + SOIRG;
   double SLRG = SWCON + SORG;
   int lent = len * (len + 1) / 2;
   int counter_1 = 0;
   int counter_2 = 0;
   arma::vec SW3(lent);
   arma::vec SG3(lent);
   arma::vec SO3(lent);
   arma::mat mat2(lent,4);
   double temp, KRW, KRG, KROW, KROG;
   for (int i = 0; i < len; i++) {
      for (int j = 0; j < (len - i); j++) {
         SW3(counter_2 + j) = SW(i);
         SG3(counter_2 + j) = SG(j);
         SO3(counter_2 + j) = 1 - SW(i) - SG(j);
         mat2(counter_2 + j,0) = SW3(counter_2 + j);
         mat2(counter_2 + j,1) = SG3(counter_2 + j);
         mat2(counter_2 + j,2) = SO3(counter_2 + j);
         if (SW3(counter_2 + j)  <= SWCRIT) {
            KRW = 0.;
         } else if (SW3(counter_2 + j)  >= (1. - SOIRW)) {
            KRW = KRWIRO;
         } else {
            KRW = KRWIRO * std::pow(((SW3(counter_2 + j)  - SWCRIT) / (1. - SWCRIT - SOIRW)), NW);
            if (std::abs(KRW) < threshold) {
               KRW = 0.0;
            }
         }
         if (SG3(counter_2 + j) <= SGCRIT) {
            KRG = 0.;
         } else if (SG3(counter_2 + j) >= (1. - SLCON)) {
            KRG = KRGCL;
         } else {
            KRG = KRGCL * std::pow(((SG3(counter_2 + j) - SGCRIT) / (1. - SGCRIT - SLCON)), NG);
            if (std::abs(KRG) < threshold) {
               KRG = 0.0;
            }
         }
         temp = 1. - SW3(counter_2 + j) ;
         if (temp <= SORW) {
            KROW = 0.;
         } else if (temp >= (1. - SWCON)) {
            KROW = KROCW;
         } else {
            KROW = KROCW * std::pow(((temp - SORW) / (1. - SWCON - SORW)), NOW);
            if (std::abs(KROW) < threshold) {
               KROW = 0.0;
            }
         }
         temp = 1. - SG3(counter_2 + j) ;
         if (temp <= SLRG) {
            KROG = 0.;
         } else if (temp >= (1. - SGCON)) {
            KROG = KROGCG;
         } else {
            KROG = KROGCG * std::pow(((temp - SLRG) / (1. - SGCON - SLRG)), NOG);
            if (std::abs(KROG) < threshold) {
               KROG = 0.0;
            }
         }
         mat2(counter_2 + j,3) = KROCW * ((KROW / KROCW + KRW) * (KROG / KROCW + KRG) - KRW - KRG);
         if (mat2(counter_2 + j,3) < 0.) {
            mat2(counter_2 + j,3) = 0.;
         }
      }
      counter_1 = len - i;
      counter_2 = counter_2 + counter_1;
   }
   NumericMatrix mat22(lent,4);
   mat22 = Rcpp::wrap(mat2);
   CharacterVector colname(4);
   colname = {"Sw", "Sg", "So", "Kro"};
   colnames(mat22) = colname;
   return(mat22);
}



//******************************************************************************
// 3-PHASE BAKER

// [[Rcpp::export]]
double krow2p_BC(double SW, double SWCON, double SORW, double KROCW, double NOW) {

   const double threshold = 1e-6;
   double SO = 1. - SW;
   double KROW;
   if (SO <= SORW) {
      KROW = 0.;
   } else if (SO >= (1 - SWCON)) {
      KROW = KROCW;
   } else {
      KROW = KROCW * std::pow(((SO - SORW) / (1 - SWCON - SORW)), NOW);
      if (std::abs(KROW) < threshold) {
         KROW = 0.0;
      }
   }
   return KROW;
}

// [[Rcpp::export]]
double krgl2p_BC(double SG, double SWCON, double SORG, double SGCON, double KROGCG, double NOG) {

   const double threshold = 1e-6;
   double SLRG = SWCON + SORG;
   double SL = 1. - SG;
   double KROG;
   if (SL <= SLRG) {
      KROG = 0;
   } else if (SL >= (1 - SGCON)) {
      KROG = KROGCG;
   } else {
      KROG = KROGCG * std::pow(((SL - SLRG) / (1 - SGCON - SLRG)), NOG);
      if (std::abs(KROG) < threshold) {
         KROG = 0.0;
      }
   }
   return KROG;
}


// //' @export
// // [[Rcpp::export]]
// arma::vec fn(arma::vec s, double SW3t, double SG3t, double SWCON, double SORW, double SORG, double SGCON, double KROCW, double KROGCG, double NOW, double NOG) {
//
//    arma::vec f(2);
//    double SW = std::pow(s(0),2.);
//    double SG = std::pow(s(1),2.);
//    f(0) = krow2p_BC(SW, SWCON, SORW, KROCW, NOW) - krgl2p_BC(SG, SWCON, SORG, SGCON, KROGCG, NOG);
//    f(1) = (SW3t - SWCON) * (SG3t - SGCON) - (SW3t - SW) *  (SG3t - SG);
//    // f(1) <- (SWCON - SW) * (SG3t - SGCON) - (SW3t - SW) *  (SG - SGCON);
//    return f;
// }


// // [[Rcpp::export]]
// arma::vec Sys_Equ(double SW3t, double SG3t, double SWCON, double SORW, double SORG, double SGCON, double KROCW, double KROGCG, double NOW, double NOG) {
//
//    Environment pkg = Environment::namespace_env("nleqslv");
//    Function func_nleqslv = pkg["nleqslv"];
//
//    arma::vec f(2);
//    arma::vec s = {std::sqrt(0.66 * (SWCON + (1 - SORW))), std::sqrt(0.66 * (SGCON + (1 - SORG)))};
//    List results(8);
//    f = fn(s, SW3t, SG3t, SWCON, SORW, SORG, SGCON, KROCW, KROGCG, NOW, NOG);
//    results = func_nleqslv(Named("x", s), Named("fn", f), Named("SW3t", SW3t), Named("SG3t", SG3t), Named("SWCON", SWCON),
//                           Named("SORW", SORW), Named("SORG", SORG), Named("SGCON", SGCON), Named("KROCW", KROCW),
//                           Named("KROGCG", KROGCG), Named("NOW", NOW), Named("NOG", NOG), Named("method", "Broyden"));
//    arma::vec res =  results(0);
//    return res;
// }


// //' Generate a matrix of three-phase relative permeability data for the water-gas-oil system using Baker's linear model
// //'
// //' The 'kr3p_Baker()' creates a table of three-phase oil relative permeability data for water, gas, and oil saturation values between zero and one. This model reads the water, and gas saturation values in the three phase region as saturation inputs into two phase relative permeability models
// //' @param SWCON connate water saturation, fraction
// //' @param SWCRIT critical water saturation, fraction
// //' @param SOIRW irreducible oil saturation, fraction
// //' @param SORW residual oil saturation, fraction
// //' @param SOIRG irreducible oil saturation, fraction
// //' @param SORG residual oil saturation, fraction
// //' @param SGCON connate gas saturation, fraction
// //' @param SGCRIT critical gas saturation, fraction
// //' @param KRWIRO water relative permeability at irreducible oil
// //' @param KROCW oil relative permeability at connate water
// //' @param KRGCL gas relative permeability at connate liquid
// //' @param KROGCG oil relative permeability at connate gas
// //' @param NW exponent term for calculating krw
// //' @param NOW exponent term for calculating krow
// //' @param NG exponent term for calculating krg
// //' @param NOG exponent term for calculating krog
// //' @param NP number of saturation points in the two-phase relative permeability tables, maximum acceptable value is 501. The number of data points in the three-phase relative permeability table is (0.5 * NP * (NP + 1))
// //' @return A matrix with water saturation, gas saturation, oil saturation, and oil relative permeability values, respectively.
// //' @examples
// //' rel_perm_wgo <- kr3p_Baker(0.15, 0.2, 0.15, 0.15, 0.2, 0.2, 0.05, 0.05, 0.4, 1, 0.3, 1, 3, 2, 4, 2.5, 101)
// //' @export
// // [[Rcpp::export]]
// arma::mat kr3p_Baker(double SWCON, double SWCRIT, double SOIRW, double SORW,
//                              double SOIRG, double SORG, double SGCON, double SGCRIT,
//                              double KRWIRO, double KROCW, double KRGCL, double KROGCG,
//                              double NW, double NOW, double NG, double NOG, int NP) {
//
//    // Environment pkg = Environment::namespace_env("nleqslv");
//    // Function func_nleqslv = pkg["nleqslv"];
//    const double threshold = 1e-6;
//    int len = NP;
//    if (len > 501) {
//       len = 501;
//    }
//    int lent = len * (len + 1) / 2;
//    int counter_1 = 0;
//    int counter_2 = 0;
//    arma::vec SW = arma::regspace(0., 1. / (len - 1), 1.);
//    arma::vec SG = arma::regspace(0., 1. / (len - 1), 1.);
//    arma::mat mat2(lent,4);
//    double SW3t, SO3t, SG3t, SOr, SW2;
//    double m = (SORG - (SORW - SGCON)) / (SWCON - (1 - SORW));
//    double b = SORG - m * SWCON;
//    arma::vec s(2);
//    List res(8);
//    for (int i = 0; i < len; i++) {
//       for (int j = 0; j < (len - i); j++) {
//          mat2(counter_2 + j,0) = SW(i);
//          mat2(counter_2 + j,1) = SG(j);
//          mat2(counter_2 + j,2) = 1 - SW(i) - SG(j);
//          SW3t = SW(i);
//          SG3t = SG(j);
//          SO3t = 1 - SW(i) - SG(j);
//          if ((SW3t > SWCON) & (SW3t < (1 - SORW)) & (SG3t < (1 - SORG - SWCON)) & (SG3t > SGCON)) {
//             SOr = m * SW3t + b;
//             if (SO3t <= SOr) {
//                mat2(counter_2 + j,3) =  0;
//             } else {
//                s = Sys_Equ(SW3t, SG3t, SWCON, SORW, SORG, SGCON, KROCW, KROGCG, NOW, NOG);
//                SW2 = std::pow(s(0),2);
//                mat2(counter_2 + j,3) = krow2p_BC(SW2, SWCON, SORW, KROCW, NOW);
//             }
//          } else if (SW3t <= SWCON) {
//             mat2(counter_2 + j,3) = krgl2p_BC(SG3t, SWCON, SORG, SGCON, KROGCG, NOG);
//          } else if (SW3t >= (1 - SORW)) {
//             mat2(counter_2 + j,3) = krow2p_BC(SW3t, SWCON, SORW, KROCW, NOW);
//          } else if (SG3t <= SGCON) {
//             mat2(counter_2 + j,3) = krow2p_BC(SW3t, SWCON, SORW, KROCW, NOW);
//          } else {
//             mat2(counter_2 + j,3) = krgl2p_BC(SG3t, SWCON, SORG, SGCON, KROGCG, NOG);
//          }
//          if (mat2(counter_2 + j,3) <= threshold) {
//             mat2(counter_2 + j,3) = 0.;
//          }
//       }
//       counter_1 = len - i;
//       counter_2 = counter_2 + counter_1;
//    }
//    return(mat2);
// }
