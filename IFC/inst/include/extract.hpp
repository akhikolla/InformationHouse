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

#ifndef IFC_EXTRACT_HPP
#define IFC_EXTRACT_HPP

#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include "utils.hpp"
#include "color.hpp"
#include "decomp.hpp"
#include "tiff.hpp"
#include "resize.hpp"
#include "matrix_logic.hpp"
using namespace Rcpp;

//' @title Matrix Normalization
//' @name cpp_normalize
//' @description
//' Normalizes a finite matrix to [0,1]
//' @param mat a finite NumericMatrix.
//' @param input_range a finite NumericVector, sets the range of the input intensity values. Default is c(0,4095).
//' values outside this range are clipped.
//' @param full_range if full_range is TRUE, then input_range will be set to 'c(0,4095)' and gamma forced to 1. Default is false.
//' @param force_range if force_range is TRUE, then input_range will be adjusted to mat range in [-4095, +inf] and gamma forced to 1. Default is false.\cr
//' Note that this parameter takes the precedence over 'input_range' and 'full_range'.
//' @param gamma correction. Default is 1, for no correction.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix hpp_normalize (const Rcpp::NumericMatrix mat, 
                                   const Rcpp::NumericVector input_range = NumericVector::create(0.0,4095.0),
                                   const bool full_range = false,
                                   const bool force_range = false, 
                                   const double gamma = 1.0) {
  Rcpp::NumericMatrix out = Rcpp::no_init_matrix(mat.nrow(), mat.ncol());
  Rcpp::NumericVector ran(2);
  double gam = gamma;
  if(force_range) {
    gam = 1.0;
    ran = hpp_check_range(hpp_check_range(mat)); // twice to ensure that mat is at least of length 2
  } else {
    if(full_range) {
      gam = 1.0;
      ran[0] = 0.0;
      ran[1] = 4095.0;
    } else {
      ran = hpp_check_range(hpp_check_range(input_range)); // twice to ensure that input_range is at least of length 2
    }
  }
  
  double diff = ran[1] - ran[0];
  if(gamma == 1.0) {
    for(R_len_t i = 0; i < mat.size(); i++) {
      if(mat[i] <= ran[0]) {
        out[i] = 0.0;
        continue;
      }
      if(mat[i] >= ran[1]) {
        out[i] = 1.0;
        continue;
      }
      out[i] = (mat[i] - ran[0])/diff;
    }
  } else {
    for(R_len_t i = 0; i < mat.size(); i++) {
      if(mat[i] <= ran[0]) {
        out[i] = 0.0;
        continue;
      }
      if(mat[i] >= ran[1]) {
        out[i] = 1.0;
        continue;
      }
      out[i] = std::pow((mat[i] - ran[0])/diff, gam);
    }
  }
  if(mat.hasAttribute("mask")) out.attr("mask") = mat.attr("mask");
  return out;
}

//' @title Matrix Cleanser
//' @name cpp_cleanse
//' @description
//' Replaces values in matrix mat according to mask msk.
//' Depending of add_noise parameter values will be replaced with noise or not.
//' @param mat a NumericMatrix.
//' @param msk a IntegerMatrix.
//' @param add_noise logical, if true adds normal noise.
//' Rcpp::Rf_rnorm(bg, sd) function is used. Default is true.
//' @param bg double, mean value of the background added if add_noise is true. Default is 0.0
//' @param sd double, standard deviation of the background added if add_noise is true. Default is 0.0
//' @return a NumericMatrix with replaced according to msk.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix hpp_cleanse (const Rcpp::NumericMatrix mat, 
                                 const Rcpp::IntegerMatrix msk,
                                 const bool add_noise = true, 
                                 const double bg = 0.0, const double sd = 0.0) {
  if(!(msk.ncol() == mat.ncol()) && (msk.nrow() == mat.nrow())) Rcpp::stop("hpp_cleanse: 'mat' and 'msk' should have same dimensions");
  Rcpp::NumericMatrix out = Rcpp::no_init_matrix(mat.nrow(), mat.ncol());
  if(add_noise) {
    for(R_len_t i = 0; i < out.size(); i++) out[i] = msk[i] ? Rf_rnorm(bg, sd):mat[i];
  } else {
    for(R_len_t i = 0; i < out.size(); i++) out[i] = msk[i] ? bg:mat[i];
  }
  out.attr("mask") = msk;
  return out;
}

//' @title Equal Sized Matrix to Matrix Writer According to Mask
//' @name cpp_mask
//' @description
//' Writes matrix B in matrix A according to mask.
//' If mask is not 0 B is written, A otherwise.
//' @param A a NumericMatrix.
//' @param B a NumericMatrix.
//' @param mask a NumericMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix hpp_mask (const Rcpp::NumericMatrix A,
                              const Rcpp::NumericMatrix B,
                              const Rcpp::NumericMatrix mask) {
  R_len_t ar = A.nrow(), ac = A.ncol();
  if((B.ncol() != ac) || (mask.ncol() != ac) || (B.nrow() != ar) || (mask.nrow() != ar) ) Rcpp::stop("hpp_mask: 'A', 'B' and 'mask' should have same dimensions");
  Rcpp::NumericMatrix out(ar, ac);
  for(R_len_t i = 0; i < A.size(); i++) out[i] = mask[i] ? B[i] : A[i];
  return out;
}

//' @title Matrix to Matrix Writer According to Mask with Offsets
//' @name cpp_mark
//' @description
//' Writes matrix B in matrix A according to mask.
//' @param A a NumericMatrix.
//' @param B a NumericMatrix.
//' @param mask a NumericMatrix.
//' @param xoff x offset in A to start writing B.
//' @param yoff x offset in A to start writing B.
//' @param invert a logical. Default is false.
//' When false, the default, values of B are written into A when mask is not 0.
//' When true, values of 1-B are written into A when mask is not 0.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix hpp_mark (const Rcpp::NumericMatrix A,
                              const Rcpp::NumericMatrix B,
                              const Rcpp::NumericMatrix mask,
                              const R_len_t xoff = 0,
                              const R_len_t yoff = 0,
                              const bool invert = false) {
  R_len_t bc = B.ncol();
  R_len_t br = B.nrow();
  if((A.ncol() < (bc + xoff)) || (A.nrow() < (br + yoff))) Rcpp::stop("hpp_mark: A should be at least of same dimensions as 'B' + 'offsets'");
  if((mask.ncol() < bc) || (mask.nrow() < br)) Rcpp::stop("hpp_mark: 'mask' should be at least of same dimensions as 'B'");
  Rcpp::NumericMatrix out = Rcpp::clone(A);
  if(invert) {
    for(R_len_t y = 0; y < br; y++) for(R_len_t x = 0; x < bc; x++) if(mask(y,x)) out(y+yoff,x+xoff) = std::fabs(1-B(y,x));
  } else {
    for(R_len_t y = 0; y < br; y++) for(R_len_t x = 0; x < bc; x++) if(mask(y,x)) out(y+yoff,x+xoff) = B(y,x);
  }
  return out;
}

//' @title Matrix Transformation
//' @name cpp_transform
//' @description
//' Function to normalize, colorize and add background to images.
//' @param mat NumericMatrix.
//' @param color NumericVector, whose members are h,s,v color.
//' This vector has to be named with 1st name being the name of this color.
//' @param msk IntegerMatrix.
//' @param size a length 2 IntegerVector, of final dimensions (height,width) of the image. Default is 0,0 for no change.
//' @param mode string, color mode export. Either "rgb", "gray" or "raw". Default is "raw".
//' @param type uint16_t image object type.
//' Rcpp::Rf_rnorm(bg, sd) function is used. Default is true.
//' @param input_range a finite NumericVector, only apply when mode is not "raw", sets the range of the input intensity values. Default is c(0,4095).
//' values exceeding this range are clipped.
//' @param add_noise bool, if true adds normal noise.
//' @param bg double, mean value of the background added if add_noise is true. Default is 0.
//' @param sd double, standard deviation of the background added if add_noise is true. Default is 0.
//' @param full_range bool, only apply when mode is not "raw", if full_range is TRUE, then 'input_range' will be set to 'c(0,4095)' and gamma forced to 1. Default is false.
//' @param force_range bool, only apply when mode is not "raw", if force_range is TRUE, then 'input_range' will be adjusted to mat range in [-4095, +inf] and gamma forced to 1. Default is false.\cr
//' Note that this parameter takes the precedence over 'input_range' and 'full_range'.
//' @param gamma correction. Default is 1, for no correction.
//' @details When add_noise is false, background will be automatically set to minimal pixel value for "masked" and "MC" removal method.\cr
//' when a mask is detected, add_noise, full_range and force_range are set to false, background mean and sd to 0, and input_range to [0,3].\cr
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::NumericVector hpp_transform(const Rcpp::NumericMatrix mat,
                                  const Rcpp::NumericVector color,
                                  const Rcpp::IntegerMatrix msk,
                                  const Rcpp::IntegerVector size = Rcpp::IntegerVector::create(0,0),
                                  const std::string mode = "raw",
                                  const uint16_t type = 2,
                                  const Rcpp::NumericVector input_range = Rcpp::NumericVector::create(0.0,4095.0),
                                  const bool add_noise = true,
                                  const double bg = 0.0,
                                  const double sd = 0.0,
                                  const bool full_range = false,
                                  const bool force_range = false,
                                  const double gamma = 1.0) {
  Rcpp::NumericMatrix foo;
  Rcpp::CharacterVector removal = Rcpp::wrap(msk.attr("removal"));
  Rcpp::CharacterVector col_name = Rcpp::wrap(color.attr("names"));
  
  if(type == 3) { // a mask is detected parameters are forced
    Rcpp::NumericVector msk_range = Rcpp::NumericVector::create(0.0, 3.0);
    if((removal[0] == "none") || (removal[0] == "raw")) {
      foo = Rcpp::clone(mat);
      foo.attr("mask") = msk;
    } else {
      foo = hpp_cleanse(mat, msk, false, 0.0, 0.0);
    } 
    foo = hpp_resize(foo, size[0], size[1], false, 0.0, 0.0);
    if(mode != "raw") {
      foo = hpp_normalize(foo, msk_range, false, false, 1.0);
    }
    if(mode == "rgb") {
      Rcpp::NumericVector bar = hpp_M_HSV2RGB(foo, color[0], color[1]);
      bar.attr("dim") = Rcpp::Dimension(foo.nrow(), foo.ncol(), 3);
      bar.attr("input_range") = msk_range;
      bar.attr("full_range") = false;
      bar.attr("force_range") = false;
      bar.attr("gamma") = 1.0;
      bar.attr("color") = Rcpp::String(col_name[0]);
      bar.attr("mode") = mode;
      bar.attr("removal") = removal;
      bar.attr("raw") = mat;
      bar.attr("BG_MEAN") = 0.0;
      bar.attr("BG_STD") = 0.0;
      bar.attr("class") = "IFC_msk";
      return bar;
    }
    foo.attr("input_range") = msk_range;
    foo.attr("full_range") = false;
    foo.attr("force_range") = false;
    foo.attr("gamma") = 1.0;
    foo.attr("color") = "Gray";
    foo.attr("mode") = mode;
    foo.attr("removal") = removal;
    foo.attr("raw") = mat;
    foo.attr("BG_MEAN") = 0.0;
    foo.attr("BG_STD") = 0.0;
    foo.attr("class") = "IFC_msk";
  } else {
    double bg_2;
    double sd_2;
    if((!add_noise) && (mode != "raw") && ((removal[0] == "masked") || (removal[0] == "MC"))) {
      bg_2 = -4096.0;
      sd_2 = 0.0;
    } else {
      bg_2 = bg;
      sd_2 = sd;
    }
    if((removal[0] == "none") || (removal[0] == "raw")) {
      foo = Rcpp::clone(mat);
      foo.attr("mask") = msk;
    } else {
      foo = hpp_cleanse(mat, msk, add_noise, bg_2, sd_2);
    } 
    foo = hpp_resize(foo, size[0], size[1], add_noise, bg_2, sd_2);
    if(mode != "raw") {
      foo = hpp_normalize(foo, input_range, full_range, force_range, gamma);
    }
    if(mode == "rgb") {
      Rcpp::NumericVector bar = hpp_M_HSV2RGB(foo, color[0], color[1]);
      bar.attr("dim") = Rcpp::Dimension(foo.nrow(), foo.ncol(), 3);
      bar.attr("input_range") = input_range;
      bar.attr("full_range") = full_range;
      bar.attr("force_range") = force_range;
      bar.attr("gamma") = gamma;
      bar.attr("color") = Rcpp::String(col_name[0]);
      bar.attr("mode") = mode;
      bar.attr("removal") = removal;
      bar.attr("raw") = mat;
      bar.attr("BG_MEAN") = bg; 
      bar.attr("BG_STD") = sd;
      bar.attr("class") = "IFC_img";
      return bar;
    }
    foo.attr("input_range") = input_range;
    foo.attr("full_range") = full_range;
    foo.attr("force_range") = force_range;
    foo.attr("gamma") = gamma;
    foo.attr("color") = "Gray";
    foo.attr("mode") = mode;
    foo.attr("removal") = removal;
    foo.attr("raw") = mat;
    foo.attr("BG_MEAN") = bg;
    foo.attr("BG_STD") = sd;
    foo.attr("class") = "IFC_img";
  }
  return foo;
}

//' @title IFC_object Extraction
//' @name cpp_extract
//' @description
//' Extracts object from ifd
//' @param fname string, path to file
//' @param ifd List, ifd information of class IFC_ifd
//' @param colors List of colors to use.
//' @param channels DataFrame, channels information
//' @param physicalChannel CharacterVector of indices for each channels
//' @param xmin NumericVector of minimal values for each channels
//' @param xmax NumericVector of maximal values for each channels
//' @param removal IntegerVector of removal method to be used for each channels
//' @param add_noise LogicalVector of whether to add_noise for each channels
//' @param full_range LogicalVector of whether to use full_range for each channels
//' @param force_range LogicalVector of whether to use force_range for each channels
//' @param gamma NumericVector of the gamma for each channels
//' @param chan_to_extract IntegerVector, channels to extract
//' @param extract_msk uint8_t, type of masked to extract.\cr
//' - 0: no mask\cr
//' - 1: at least one raw mask\cr
//' - 2: at least one clipped\cr
//' - 3: at least one masked\cr
//' - 4: at least one MC
//' @param mode string, color mode export. Either "rgb", "gray" or "raw". Default is "raw".
//' @param size a length 2 IntegerVector of final dimensions (height,width) of the image. Default is 0,0 for no change.\cr
//' @param verbose bool, whether to display information (use for debugging purpose). Default is false.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::List hpp_extract (const std::string fname,
                        const Rcpp::List ifd,
                        const Rcpp::List colors,
                        const Rcpp::CharacterVector physicalChannel,
                        const Rcpp::NumericVector xmin,
                        const Rcpp::NumericVector xmax,
                        const Rcpp::IntegerVector removal,
                        const Rcpp::LogicalVector add_noise,
                        const Rcpp::LogicalVector full_range,
                        const Rcpp::LogicalVector force_range,
                        const Rcpp::NumericVector gamma,
                        const Rcpp::IntegerVector chan_to_extract,
                        const uint8_t extract_msk = 0,
                        const std::string mode = "raw",
                        const Rcpp::IntegerVector size = Rcpp::IntegerVector::create(0,0),
                        const bool verbose = false) {
  R_len_t nb_channels = physicalChannel.size();
  Rcpp::List infos = ifd["infos"];
  R_len_t iml = infos["IMAGE_LENGTH"];
  R_len_t imw = infos["IMAGE_WIDTH"];
  uint16_t typ = infos["TYPE"];
  uint32_t off = infos["STRIP_OFFSETS"];
  uint32_t byt = infos["STRIP_BYTE_COUNTS"];
  uint32_t com = infos["COMPRESSION"];
  Rcpp::NumericVector bg_2;
  Rcpp::NumericVector sd_2;
  Rcpp::List out(chan_to_extract.size());
  
  switch(typ) {
  case 1: { // 1st offset is detected, nothing to extract
    for(R_len_t i = 0; i < chan_to_extract.size(); i++) {
    out[i] = Rcpp::NumericMatrix(iml, imw);
  }
    return out;
    break;
  }
  case 2: { // an image is detected use background mean and sd
    bg_2 = Rcpp::clone(Rcpp::as<Rcpp::NumericVector>(infos["BG_MEAN"]));
    sd_2 = Rcpp::clone(Rcpp::as<Rcpp::NumericVector>(infos["BG_STD"]));
    break;
  }
  case 3: { // a mask is detected some parameters are forced
    bg_2 = Rcpp::rep(0.0, nb_channels);
    sd_2 = Rcpp::rep(0.0, nb_channels);
    break;
  }
  default: { // not allowed type
    Rcpp::stop("hpp_extract: trying to extract a unknow object");
  }
  }
  
  // extract image
  Rcpp::List img = hpp_decomp(fname, off, byt, 
                        imw, iml, nb_channels,
                        0, com, verbose);
  Rcpp::IntegerMatrix msk_init(iml, imw / nb_channels);
  msk_init.fill(0);
  Rcpp::IntegerMatrix msk = Rcpp::clone(msk_init);
  Rcpp::IntegerMatrix MC = Rcpp::clone(msk_init);
  
  if(extract_msk > 0) {
    Rcpp::List masks;
    if(typ == 2) {
      Rcpp::List msk_ifd = hpp_getTAGS(fname, ifd["next_IFD_offset"], verbose, 8, true)["infos"];
      masks = hpp_decomp(fname, msk_ifd["STRIP_OFFSETS"], msk_ifd["STRIP_BYTE_COUNTS"], 
                         imw, iml, nb_channels, 
                         0, msk_ifd["COMPRESSION"], verbose);
      
    } else {
      masks = clone(img);
    }
    if(extract_msk == 4) {
      for(R_len_t i = 0; i < masks.length(); i++) {
        Rcpp::IntegerMatrix CUR_M = Rcpp::clone(Rcpp::as<Rcpp::IntegerMatrix>(masks[i]));
        for(R_len_t i_row = 0; i_row < iml; i_row ++) {
          MC(i_row, Rcpp::_) = (CUR_M(i_row, Rcpp::_) == 1) | (MC(i_row, Rcpp::_) == 1);
        }
      }
      for(R_len_t i_row = 0; i_row < iml; i_row ++) {
        MC(i_row, Rcpp::_) = (MC(i_row, Rcpp::_) == 0);
      }
    }
    masks.attr("names") = physicalChannel;
    
    // transform extracted image according to user's settings with mask removal
    for(R_len_t i = 0; i < chan_to_extract.size(); i++) {
      msk = Rcpp::clone(msk_init); // mandatory to be sure that it is new matrix each time
      uint8_t chan_idx = chan_to_extract[i];
      std::string cur_chan = as<std::string>(physicalChannel[chan_idx]);
      switch(removal[chan_idx]) {
      case 0: {
        msk.attr("removal") = "none";
        break;
      }
      case 1: {
        msk = Rcpp::clone(Rcpp::as<Rcpp::IntegerMatrix>(masks[cur_chan]));
        msk.attr("removal") = "raw";
        break;
      }
      case 2: {
        Rcpp::IntegerMatrix CUR_M = Rcpp::clone(Rcpp::as<Rcpp::IntegerMatrix>(masks[cur_chan]));
        for(R_len_t i_row = 0; i_row < iml; i_row ++) {
          msk(i_row, Rcpp::_) = CUR_M(i_row, Rcpp::_) > 1;
        }
        msk.attr("removal") = "clipped";
        break;
      }
      case 3: {
        Rcpp::IntegerMatrix CUR_M = Rcpp::clone(Rcpp::as<Rcpp::IntegerMatrix>(masks[cur_chan]));
        for(R_len_t i_row = 0; i_row < iml; i_row ++) {
          msk(i_row, Rcpp::_) = CUR_M(i_row, Rcpp::_) != 1;
        }
        msk.attr("removal") = "masked";
        break;
      }
      case 4: {
        msk = MC;
        msk.attr("removal") = "MC";
        break;
      }
      }
      
      out[i] = hpp_transform(img[chan_idx],
                             colors[chan_idx],
                                   msk,
                                   size,
                                   mode,
                                   typ,
                                   NumericVector::create(xmin[chan_idx],xmax[chan_idx]),
                                   add_noise[chan_idx],
                                            bg_2[chan_idx],
                                                sd_2[chan_idx],
                                                    full_range[chan_idx],
                                                              force_range[chan_idx],
                                                                         gamma[chan_idx]);
    }
  } else {
    msk.attr("removal") = "none";
    for(R_len_t i = 0; i < chan_to_extract.size(); i++) {  
      // transform extracted image according to user's settings without mask removal
      uint8_t chan_idx = chan_to_extract[i];
      out[i] = hpp_transform(img[chan_idx],
                             colors[chan_idx],
                                   msk,
                                   size,
                                   mode,
                                   typ,
                                   NumericVector::create(xmin[chan_idx],xmax[chan_idx]),
                                   add_noise[chan_idx],
                                            bg_2[chan_idx],
                                                sd_2[chan_idx],
                                                    full_range[chan_idx],
                                                              force_range[chan_idx],
                                                                         gamma[chan_idx]);
    }
  }
  return out;
}

#endif
