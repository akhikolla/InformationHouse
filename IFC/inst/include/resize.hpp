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

#ifndef IFC_RESIZE_HPP
#define IFC_RESIZE_HPP

#include <Rcpp.h>
using namespace Rcpp;

//' @title Header for Matrix Cropping
//' @name cpp_crop
//' @description
//' Crops mat according to new_height and new_width parameters.
//' @param mat a numeric matrix.
//' @param new_height an unsigned integer, giving the new height of returned mat. Default is 0 for no change.
//' @param new_width an unsigned integer, giving the new width of returned mat. Default is 0 for no change.
//' @return a cropped matrix.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix hpp_crop (Rcpp::NumericMatrix mat,
                              const R_len_t new_height = 0,
                              const R_len_t new_width = 0) {
  R_len_t img_c, img_r;
  img_c = mat.ncol();
  img_r = mat.nrow();
  // new dimensions are larger than original dimensions, no need to crop
  if((img_c <= new_width) && (img_r <= new_height)) return mat;
  // new width is larger than original width and new_height is 0, no need to crop
  if((img_c <= new_width) && (new_height == 0)) return mat;
  // new height is larger than original height and new_width is 0, no need to crop
  if((img_r <= new_height) && (new_width == 0)) return mat;
  
  R_len_t ori_c, ori_r, crop_c, crop_r;
  // height resizing
  if((new_height > 0) && (new_height < img_r)) {
    // compute same amount of rows to remove on top-bottom
    ori_r = (img_r - new_height) >> 1;
    crop_r = new_height;
  } else { // no height resizing
    ori_r = 0;
    crop_r = img_r;
  }
  // width resizing
  if((new_width > 0) && (new_width < img_c)) {
    // compute same amount of cols to remove on right-left
    ori_c = (img_c - new_width) >> 1;
    crop_c = new_width;
  } else { // no width resizing
    ori_c = 0;
    crop_c = img_c;
  }
  return mat( Rcpp::Range(ori_r,crop_r + ori_r - 1) , Rcpp::Range(ori_c,crop_c + ori_c - 1) );
}

// function to expand matrix without adding noise
Rcpp::NumericMatrix hpp_expand_no_noise (const Rcpp::NumericMatrix mat,
                                         const R_len_t new_height = 0,
                                         const R_len_t new_width = 0,
                                         const double bg = 0.0) {
  R_len_t img_r, img_c;
  img_r = mat.nrow();
  img_c = mat.ncol();
  // new dimensions are smaller than original dimensions, no need to expand
  if((img_r >= new_height) && (img_c >= new_width)) return mat;
  
  R_len_t i_row, i_col, ori_r, ori_c, i_out, fin_height, fin_width;
  // final height: either original or desired height if larger
  fin_height = img_r > new_height ? img_r : new_height;
  // final width: either original or desired width if larger
  fin_width = img_c > new_width ? img_c : new_width;
  // compute padding
  ori_r = (fin_height-img_r) >> 1;
  ori_c = (fin_width-img_c) >> 1;
  
  // create output matrix with expanded dimension and fill with bg
  Rcpp::NumericMatrix out(fin_height, fin_width);
  out.fill(bg);
  
  // write mat into center of output matrix
  for(i_col = 0; i_col < img_c; i_col++) {
    i_out = (i_col + ori_c) * fin_height + ori_r;
    for(i_row = 0; i_row < img_r; i_row++) {
      out[i_out++] = mat(i_row, i_col);
    }
  }
  return out;
}

// function to expand matrix with new rows adding noisy padding bg
Rcpp::NumericMatrix hpp_expand_row (const Rcpp::NumericMatrix mat,
                                    const R_len_t new_height = 0,
                                    const double bg = 0.0, const double sd = 0.0) {
  R_len_t img_r = mat.nrow();
  // new dimensions are smaller than original dimensions, no need to expand
  if(img_r >= new_height) return mat;
  
  R_len_t img_c, ori_r, i;
  img_c = mat.ncol();
  // compute row padding
  ori_r = (new_height-img_r) >> 1;
  
  // create output matrix with expanded rows
  Rcpp::NumericMatrix out(new_height, img_c);
  
  // write top padding
  for(i = 0; i < ori_r; i++) out(i, Rcpp::_) = Rcpp::rnorm(img_c, bg, sd);
  // write mat into center of output matrix
  for(; i < (ori_r + img_r); i++) out(i, Rcpp::_) = mat(i - ori_r, Rcpp::_);
  // write bottom padding
  for(; i < new_height; i++) out(i, Rcpp::_) = Rcpp::rnorm(img_c, bg, sd); 
  return out;
}

// function to expand matrix with new columns adding noisy padding bg
Rcpp::NumericMatrix hpp_expand_col (const Rcpp::NumericMatrix mat,
                                    const R_len_t new_width = 0,
                                    const double bg = 0.0, const double sd = 0.0) {
  R_len_t img_c = mat.ncol();
  // new dimensions are smaller than original dimensions, no need to expand
  if(img_c >= new_width) return mat;
  
  R_len_t img_r, ori_c, i;
  img_r = mat.nrow();
  // compute row padding
  ori_c = (new_width-img_c) >> 1;
  
  // create output matrix with expanded rows
  Rcpp::NumericMatrix out(img_r, new_width);
  
  // write left padding
  for(i = 0; i < ori_c; i++) out(Rcpp::_, i) = Rcpp::rnorm(img_r, bg, sd);
  // write mat into center of output matrix
  for(; i < (ori_c + img_c); i++) out(Rcpp::_, i) = mat(Rcpp::_, i - ori_c);
  // write right padding
  for(; i < new_width; i++) out(Rcpp::_, i) = Rcpp::rnorm(img_r, bg, sd);
  
  return out;
}

// function to expand matrix with padding noisy bg
Rcpp::NumericMatrix hpp_expand_w_noise (const Rcpp::NumericMatrix mat, 
                                        const R_len_t new_height = 0,
                                        const R_len_t new_width = 0,
                                        const double bg = 0.0, const double sd = 0.0) {
  Rcpp::NumericMatrix M0 = hpp_expand_col(mat, new_width, bg, sd);
  return hpp_expand_row(M0, new_height, bg, sd);
}

//' @title Header for Matrix Resizing
//' @name cpp_resize
//' @description
//' Resizes mat according to new_height and new_width parameters.
//' @param mat a numeric matrix.
//' @param new_height an unsigned integer, giving the new height of returned mat. Default is 0 for no change.
//' @param new_width an unsigned integer, giving the new width of returned mat. Default is 0 for no change.
//' @param add_noise logical, if true adds normal noise when at least one new dimension is larger than original mat dimensions 
//' Rcpp::rnorm() function is used. Default is true.
//' @param bg double, mean value of the background added if add_noise is true. Default is 0.
//' @param sd double, standard deviation of the background added if add_noise is true. Default is 0.
//' @return a resized matrix with padding background if new_height or new_width is larger than original mat dimensions.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix hpp_resize (const Rcpp::NumericMatrix mat, 
                                const R_len_t new_height = 0, 
                                const R_len_t new_width = 0,
                                const bool add_noise = true, 
                                const double bg = 0.0, const double sd = 0.0) {
  Rcpp::NumericMatrix crop = hpp_crop(mat, new_height, new_width);
  Rcpp::NumericMatrix out;
  if(add_noise) {
    out = hpp_expand_w_noise(crop, new_height, new_width, bg, sd);
  } else {
    out = hpp_expand_no_noise(crop, new_height, new_width, bg);
  }
  if(mat.hasAttribute("mask")) out.attr("mask") = mat.attr("mask");
  return out;
}

#endif
