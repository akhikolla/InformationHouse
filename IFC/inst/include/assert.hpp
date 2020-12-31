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

#ifndef IFC_ASSERT_HPP
#define IFC_ASSERT_HPP

#include <Rcpp.h>
using namespace Rcpp;

// combines 2 character vector into one
Rcpp::CharacterVector hpp_combine(const Rcpp::CharacterVector x,
                                  const Rcpp::CharacterVector y) {
  Rcpp::CharacterVector out(x.size() + y.size());
  R_len_t i = 0, j = 0;
  for(; i<x.size(); i++) out[i] = x[i];
  for(; j<y.size(); i++, j++) out[i] = y[j];
  return out;
}

// extracts Rtype of object
Rcpp::CharacterVector hpp_type(const RObject x) {
  switch( TYPEOF(x) ){
  case NILSXP: return "NULL"; // 0
  case SYMSXP: return "symbol"; // 1
  case LISTSXP: return "pairlist"; // 2
  case CLOSXP: return "closure"; // 3
  case ENVSXP: return "environment"; // 4
  case PROMSXP: return "promise"; // 5
  case LANGSXP: return "language"; // 6
  case SPECIALSXP: return "special"; // 7
  case BUILTINSXP: return "builtin"; // 8 
  case CHARSXP: return "char"; // 9 .....
  case LGLSXP: return "logical"; // 10
  case INTSXP: return "integer"; // 13
  case REALSXP: return "numeric"; // 14
  case CPLXSXP: return "complex"; // 15
  case STRSXP: return "character"; // 16
  case DOTSXP: return "..."; // 17
  case ANYSXP: return "any"; // 18
  case VECSXP: return "list"; // 19
  case EXPRSXP: return "expression"; // 20
  case BCODESXP: return "bytecode"; // 21
  case EXTPTRSXP: return "externalptr"; // 22 
  case WEAKREFSXP: return "weakref"; // 23
  case RAWSXP: return "raw"; // 24
  case S4SXP: return "s4"; // 25
  default: Rcpp::stop("hpp_type: not supported type in 'x'"); // no support for obj
  }
}

// checks if members of object x only contains allowed elements from 'alw';
// ONLY works with vectors: LogicalVector, IntegerVector, NumericVector, ComplexVector and CharacterVector.
bool hpp_allowed (const RObject x, const RObject alw) {
  bool out = true;
  R_len_t i, j;
  switch( TYPEOF(x) ) {
  case LGLSXP: {
    Rcpp::LogicalVector x_copy = Rcpp::as<Rcpp::LogicalVector>(x); // 10
    Rcpp::LogicalVector alw_copy = Rcpp::as<Rcpp::LogicalVector>(alw);
    Rcpp::LogicalVector b(alw_copy.size());
    for(i = 0; i < x_copy.size() && out ; i++) {
      for(j = 0; j < alw_copy.size(); j++) {
        b[j] = (alw_copy[j] == x_copy[i]);
      }
      out = any(b).is_true();
    }
  }
    break;
  case INTSXP: {
    Rcpp::IntegerVector x_copy = Rcpp::as<Rcpp::IntegerVector>(x); // 13
    Rcpp::IntegerVector alw_copy = Rcpp::as<Rcpp::IntegerVector>(alw);
    Rcpp::LogicalVector b(alw_copy.size());
    for(i = 0; i < x_copy.size() && out ; i++) {
      for(j = 0; j < alw_copy.size(); j++) {
        b[j] = (alw_copy[j] == x_copy[i]);
      }
      out = any(b).is_true();
    }
  }
    break;
  case REALSXP: {
    Rcpp::NumericVector x_copy = Rcpp::as<Rcpp::NumericVector>(x); // 14
    Rcpp::NumericVector alw_copy = Rcpp::as<Rcpp::NumericVector>(alw);
    Rcpp::LogicalVector b(alw_copy.size());
    for(i = 0; i < x_copy.size() && out ; i++) {
      for(j = 0; j < alw_copy.size(); j++) {
        b[j] = (alw_copy[j] == x_copy[i]);
      }
      out = any(b).is_true();
    }
  }
    break;
  case CPLXSXP: {
    Rcpp::ComplexVector x_copy = Rcpp::as<Rcpp::ComplexVector>(x); // 15
    Rcpp::ComplexVector alw_copy = Rcpp::as<Rcpp::ComplexVector>(alw);
    Rcpp::LogicalVector b(alw_copy.size());
    for(i = 0; i < x_copy.size() && out ; i++) {
      for(j = 0; j < alw_copy.size(); j++) {
        b[j] = (alw_copy[j] == x_copy[i]);
      }
      out = any(b).is_true();
    }
  }
    break;
  case STRSXP: {
    Rcpp::CharacterVector x_copy = Rcpp::as<Rcpp::CharacterVector>(x); // 16
    Rcpp::CharacterVector alw_copy = Rcpp::as<Rcpp::CharacterVector>(alw);
    Rcpp::LogicalVector b(alw_copy.size());
    for(i = 0; i < x_copy.size() && out ; i++) {
      for(j = 0; j < alw_copy.size(); j++) {
        b[j] = (alw_copy[j] == x_copy[i]);
      }
      out = any(b).is_true();
    }
  }
    break;
  case RAWSXP: {
    Rcpp::RawVector x_copy = Rcpp::as<Rcpp::RawVector>(x); // 24
    Rcpp::RawVector alw_copy = Rcpp::as<Rcpp::RawVector>(alw);
    Rcpp::LogicalVector b(alw_copy.size());
    for(i = 0; i < x_copy.size() && out ; i++) {
      for(j = 0; j < alw_copy.size(); j++) {
        b[j] = (alw_copy[j] == x_copy[i]);
      }
      out = any(b).is_true();
    }
  }
    break;
  default: Rcpp::stop("hpp_allowed: not supported type in 'x'"); // no support for obj
  }
  return out;
}

//' @title Input Parameters Assertive Tool
//' @name cpp_assert
//' @description
//' Ensures that x respects several parameters
//' @param len IntegerVector, of allowed length for x. Default is NULL, for not checking this parameter.
//' @param cla CharacterVector, of allowed classes of x. Default is NULL, for not checking this parameter.
//' @param typ CharacterVector, of allowed types of x. Default is NULL, for not checking this parameter.
//' @param Robject of allowed values for x (will be passed to cpp_allowed). Default is NULL, for not checking this parameter.
//' @param fun CharacterVector, function to execute when mandatory parameters are not met. Default is "stop". Allowed are "stop","warning","message","return".
//' fun is placed in cpp_assert() in order to check it is correct before being used in assert() function.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::LogicalVector hpp_assert(const RObject x,
                               const Rcpp::Nullable<Rcpp::IntegerVector> len = R_NilValue,
                               const Rcpp::Nullable<Rcpp::CharacterVector> cla = R_NilValue,
                               const Rcpp::Nullable<Rcpp::CharacterVector> typ = R_NilValue,
                               const RObject alw = R_NilValue,
                               const Rcpp::CharacterVector fun = "stop") {
  Rcpp::LogicalVector ele = Rcpp::LogicalVector::create(len.isNotNull(), cla.isNotNull(), typ.isNotNull(), !alw.isNULL());
  if(fun.size() != 1) Rcpp::stop("'fun' should be of length 1");
  if(!hpp_allowed(fun, Rcpp::CharacterVector::create("stop", "warning", "message", "return"))) Rcpp::stop("'fun' has to be either 'stop', 'warning', 'message', 'return'");
  if(any(ele).is_true()) {
    if(ele[0]) {
      Rcpp::IntegerVector len_v = Rcpp::as<Rcpp::IntegerVector>(len);
      for(R_len_t i = 0; i<len_v.size(); i++) {
        if(len_v[i] == Rf_xlength(x)) {
          ele[0] = false;
          break;
        }
      }
      
    } 
    if(ele[1]) {
      SEXP cl = x.attr("class");
      Rcpp::CharacterVector K;
      if(!Rf_isNull(cl) && (Rf_length(cl) > 0)) {
        K = hpp_combine(Rcpp::as<Rcpp::CharacterVector>(cl), hpp_type(x));
      } else {
        K = hpp_type(x);
      }
      Rcpp::CharacterVector cla_v = Rcpp::as<Rcpp::CharacterVector>(cla);
      R_len_t cla_i = 0; 
      for(R_len_t i = 0; i<cla_v.size(); i++) {
        for(R_len_t k = 0; k<K.size(); k++) {
          if(cla_v[i] == K[k]) {
            cla_i++;
            break;
          }
        }
      }
      if(cla_v.size() == cla_i) {
        ele[1] = false;
      }
    }
    if(ele[2]) {
      if(hpp_allowed(hpp_type(x), Rcpp::as<Rcpp::CharacterVector>(typ))) {
        ele[2] = false;
      }
    }
    if(ele[3]) {
      if(hpp_allowed(x, alw)) {
        ele[3] = false;
      } 
    }
  }
  return ele;
}

#endif
