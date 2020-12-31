/*
  This file is released under the GNU General Public License, Version 3, GPL-3  
  Copyright (C) 2020 Yohann Demont                                              
                                                                                
  It is part of IFC package, please cite:                                       
  -IFC: An R Package for Imaging Flow Cytometry                                 
  -YEAR: 2020                                                                   
  -COPYRIGHT HOLDERS: Yohann Demont, Gautier Stoll, Guido Kroemer,              
                      Jean-Pierre Marolleau, Loï?c Gaççon,                       
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

#ifndef IFC_MATRIX_LOGIC_HPP
#define IFC_MATRIX_LOGIC_HPP

#include <Rcpp.h>
using namespace Rcpp;

//' @title Matrix List And Logic
//' @name cpp_AND_M
//' @description
//' This function takes a list of matrices and returns the AND operation applied on these matrices.
//' @param list a list of logical matrices.
//' @return a logical matrix.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::LogicalMatrix hpp_AND_M(const Rcpp::List list) {
  R_len_t L = list.length();
  if(L < 1) Rcpp::stop("hpp_AND_M: 'list' should contain at least 1 matrix");
  Rcpp::LogicalMatrix MAT = Rcpp::clone(Rcpp::as<Rcpp::LogicalMatrix>(list[0]));
  R_len_t mat_r = MAT.nrow(), mat_c = MAT.ncol();
  if(L > 1) for(R_len_t i = 1; i < L; i++) {
    Rcpp::LogicalMatrix CUR_M = Rcpp::clone(Rcpp::as<Rcpp::LogicalMatrix>(list[i]));
    if(mat_r != CUR_M.nrow()) Rcpp::stop("hpp_AND_M: 'All matrices in 'list' should have same number of rows/columns");
    if(mat_c != CUR_M.ncol()) Rcpp::stop("hpp_AND_M: 'All matrices in 'list' should have same number of rows/columns");
    for(R_len_t i_col = 0; i_col < mat_c; i_col++) {
      MAT(Rcpp::_, i_col) = CUR_M(Rcpp::_, i_col) & MAT(Rcpp::_, i_col);
    }
  }
  return MAT;
}

//' @title Matrix List Or Logic
//' @name cpp_OR_M
//' @description
//' This function takes a list of matrices and returns the OR operation applied on these matrices.
//' @param list a list of logical matrices.
//' @return a logical matrix.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::LogicalMatrix hpp_OR_M(const Rcpp::List list) {
  R_len_t L = list.length();
  if(L < 1) Rcpp::stop("hpp_OR_M: 'list' should contain at least 1 matrix");
  Rcpp::LogicalMatrix MAT = Rcpp::clone(Rcpp::as<Rcpp::LogicalMatrix>(list[0]));
  R_len_t mat_r = MAT.nrow(), mat_c = MAT.ncol();
  if(L > 1) for(R_len_t i = 1; i < L; i++) {
    Rcpp::LogicalMatrix CUR_M = Rcpp::clone(Rcpp::as<Rcpp::LogicalMatrix>(list[i]));
    if(mat_r != CUR_M.nrow()) Rcpp::stop("hpp_OR_M: 'All matrices in 'list' should have same number of rows/columns");
    if(mat_c != CUR_M.ncol()) Rcpp::stop("hpp_OR_M: 'All matrices in 'list' should have same number of rows/columns");
    for(R_len_t i_col = 0; i_col < mat_c; i_col++) {
      MAT(Rcpp::_, i_col) = CUR_M(Rcpp::_, i_col) | MAT(Rcpp::_, i_col);
    }
  }
  return MAT;
}

//' @title Matrix Neg Logic
//' @name cpp_NEG_M
//' @description
//' This function takes a logical matrix and returns its negation.
//' @param mat LogicalMatrix.
//' @return a logical matrix.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::LogicalMatrix hpp_NEG_M(const Rcpp::LogicalMatrix mat) {
  Rcpp::LogicalMatrix OUT_M = Rcpp::no_init_matrix(mat.nrow(), mat.ncol());
  for(R_len_t i_col = 0; i_col < mat.ncol(); i_col++) {
    OUT_M(Rcpp::_, i_col) = !OUT_M(Rcpp::_, i_col);
  }
  return OUT_M;
}

#endif
