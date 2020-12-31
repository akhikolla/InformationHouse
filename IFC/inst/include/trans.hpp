/*
  This file is released under the GNU General Public License, Version 3, GPL-3  
  Copyright (C) 2020 Yohann Demont                                              
                                                                                
  It is part of IFC package, please cite:                                       
  -IFC: An R Package for Imaging Flow Cytometry                                 
  -YEAR: 2020                                                                   
  -COPYRIGHT HOLDERS: Yohann Demont, Gautier Stoll, Guido Kroemer,              
                      Jean-Pierre Marolleau, Loïc Garçon,                       
                      INSERM, UPD, CHU Amiens                                   
                                                                                
                                                                                
  DISCLAIMER:                                                                   
  -You are using this package on your own risk!                                 
  -We do not guarantee privacy nor confidentiality.                             
  -This program is distributed in the hope that it will be useful, but WITHOUT  
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or         
  FITNESS FOR A PARTICULAR PURPOSE. In no event shall the copyright holders or  
  contributors be liable for any direct, indirect, incidental, special,         
  exemplary, or consequential damages (including, but not limited to,           
  procurement of substitute goods or services; loss of use, data, or profits;   
  or business interruption) however caused and on any theory of liability,      
  whether in contract, strict liability, or tort (including negligence or       
  otherwise) arising in any way out of the use of this software, even if        
  advised of the possibility of such damage.                                    
                                                                                
  You should have received a copy of the GNU General Public License             
  along with IFC. If not, see <http://www.gnu.org/licenses/>.                   
*/

#ifndef IFC_TRANS_HPP
#define IFC_TRANS_HPP

#include <Rcpp.h>
#include "utils.hpp"
using namespace Rcpp;

//' @title Smooth LinLog Transformation with Rcpp
//' @name cpp_smoothLinLog
//' @description
//' Takes a numeric vector and return its transformation:
//' - to linear, if abs(x) < hyper.
//' - to log, if abs(x) > hyper.
//' @param x NumericVector.
//' @param hyper double, value where transition between Lin/Log is applied.
//' @param base double, base of Log scale.
//' @param lin_comp double, value that is used to smooth transition between Lin/Log.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::NumericVector hpp_smoothLinLog (const Rcpp::NumericVector x,
                                      const double hyper = 1000.0,
                                      const double base = 10.0,
                                      const double lin_comp = 2.302585) {
  if(nNotisNULL(x)) {
    R_len_t L = x.size();
    double K = std::log(base) / lin_comp;
    Rcpp::NumericVector B = K + K * Rcpp::log((Rcpp::abs(x))/hyper);
    Rcpp::NumericVector xx = Rcpp::no_init_vector(L);
    for(R_len_t i = 0;i < L; i++) {
      if( std::fabs(x(i)) <= hyper ) {
        xx(i) = x(i) * K / hyper;
      } else {
        xx(i) = B(i);
        if(x(i)<0) xx(i) *= -1;
      }
    }
    return xx;
  } else {
    return x;
  }
}

//' @title Inverse Smooth LinLog Transformation with Rcpp
//' @name cpp_inv_smoothLinLog
//' @description
//' Takes a numeric vector and return its transformation:
//' - to linear, if abs(x) < log(base) / lin_comp.
//' - to exp, if abs(x) > log(base) / lin_comp.
//' @param x NumericVector.
//' @param hyper double, value where transition between Lin/Log is applied.
//' @param base double, base of Log scale.
//' @param lin_comp double, value that is used to smooth transition between Lin/Log.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::NumericVector hpp_inv_smoothLinLog (const Rcpp::NumericVector x,
                                          const double hyper = 1000.0, 
                                          const double base = 10.0, 
                                          const double lin_comp = 2.302585) {
  if(nNotisNULL(x)) {
    R_len_t L = x.size();
    double K = std::log(base) / lin_comp;
    Rcpp::NumericVector B = hyper * Rcpp::exp((Rcpp::abs(x))/K-1);
    Rcpp::NumericVector xx = Rcpp::no_init_vector(L);
    for(R_len_t i = 0;i < L; i++) {
      if( std::fabs(x(i)) > K ) {
        xx(i) = B(i);
        if(x(i)<0) xx(i) *= -1;
      } else {
        xx(i) = x(i) * hyper / K;
      }
    }
    return xx;
  } else {
    return x;
  }
}

//' @title Int32 to Uint32 32bits Conversion
//' @name cpp_int32_to_uint32
//' @description
//' Converts 32bits integer from signed to unsigned
//' @param x int32_t.
//' @keywords internal
////' @export
// [[Rcpp::export]]
uint32_t hpp_int32_to_uint32 (const int32_t x) {
  uint32_t out = x;
  return out;
}

//' @title Uint32 to Int32 32bits Conversion
//' @name cpp_uint32_to_int32
//' @description
//' Converts 32bits integer from unsigned to signed
//' @param x uint32_t.
//' @keywords internal
////' @export
// [[Rcpp::export]]
int32_t hpp_uint32_to_int32 (const uint32_t x) {
  int32_t out = x;
  return out;
}

//' @title Numeric to String Conversion
//' @name cpp_num_to_string
//' @description
//' Formats numeric to string used for features, images, ... values conversion when exporting to xml.
//' @param x a numeric vector.
//' @param precision number of significant decimal digits to keep when abs(x) < 1. Default is 15.
//' @return a string vector.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::StringVector hpp_num_to_string(const Rcpp::NumericVector x, const unsigned char precision = 16) {
  Rcpp::StringVector out(x.size());
  std::stringstream stream;
  std::string s;
  unsigned char relative;
  for(R_len_t i =0; i < x.size(); i++) {
    if(x(i) == R_PosInf) {
      out(i) = "+Inf";
      continue;
    }
    if(x(i) == R_NegInf) {
      out(i) = "-Inf";
      continue;
    }
    if(std::isnan(x[i])) { // for NaN or NA
      out(i) = "NaN";
      continue;
    }
    if(x[i] < 1) {
      relative = 0;
    } else {
      relative = std::log10(x[i]);
      if(x[i] > std::pow(10, relative)) {
        relative = relative + 1;
      } else {
        relative = relative + 2;
      }
    }
    stream << std::fixed << std::setprecision(precision - relative) << x(i);
    s = stream.str();
    s.erase(s.find_last_not_of('0') + 1, s.size());
    s.erase(s.find_last_not_of('.') + 1, s.size());
    out(i) = s;
    stream.str(std::string()); // clear stream
  }
  return out;
}

#endif
