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

#ifndef IFC_DECOMP_HPP
#define IFC_DECOMP_HPP

#include <Rcpp.h>
#include <iostream>
#include <fstream>
using namespace Rcpp;

//' @title RLE Decompression
//' @name cpp_rle_Decomp
//' @description
//' Operates RLE decompression of compressed image stored in TIFF file.
//' @param fname string, path to file.
//' @param offset uint32_t, position of the beginning of compressed image.
//' @param nbytes uint32_t, number of bytes of compressed image.
//' @param imgWidth R_len_t, Width of the decompressed image. Default is 1.
//' @param imgHeight R_len_t, Height of the decompressed image. Default is 1.
//' @param nb_channels R_len_t, number of channels of the decompressed image. Default is 1.
//' @param removal uint8_t, object removal method. Default is no removal. Otherwise, if\cr
//' -1, for clipped removal: height OR width clipped pixels will be set to -1.\cr
//' -2, height clipped removal: height clipped pixels will be set to -1.\cr
//' -3, width clipped removal: width clipped pixels will be set to -1.\cr
//' -4, only keep background: background pixels will be set to 1 and all others to -1.\cr
//' -5, only keep foreground: foreground pixels will be set to 1 and all others to -1
//' @param verbose bool, whether to display information (use for debugging purpose). Default is false.
//' @details
//' BSD implementations of Bio-Formats readers and writers
//' %%
//' Copyright (C) 2005 - 2017 Open Microscopy Environment:
//'   - Board of Regents of the University of Wisconsin-Madison
//'   - Glencoe Software, Inc.
//'   - University of Dundee
//' %%
//' Redistribution and use in source and binary forms, with or without
//' modification, are permitted provided that the following conditions are met:
//' 
//' 1. Redistributions of source code must retain the above copyright notice,
//'    this list of conditions and the following disclaimer.
//' 2. Redistributions in binary form must reproduce the above copyright notice,
//'     this list of conditions and the following disclaimer in the documentation
//'     and/or other materials provided with the distribution.
//'  
//' THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
//' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
//' IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
//' ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
//' LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
//' CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
//' SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
//' INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
//' CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
//' ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
//' POSSIBILITY OF SUCH DAMAGE.
//' @source For image decompression, Lee Kamentsky's code porting from \url{https://github.com/ome/bioformats/blob/4146b9a1797501f0fec7d6cfe69124959bff96ee/components/formats-bsd/src/loci/formats/in/FlowSightReader.java}\cr
//' cited in \url{https://linkinghub.elsevier.com/retrieve/pii/S1046-2023(16)30291-2}
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::List hpp_rle_Decomp (const std::string fname, 
                           const uint32_t offset,
                           const uint32_t nbytes,
                           const R_len_t imgWidth = 1,
                           const R_len_t imgHeight = 1,
                           const R_len_t nb_channels = 1,
                           const uint8_t removal = 0,
                           const bool verbose = false) {
  if((nb_channels * imgWidth * imgHeight) != 0) {
    Rcpp::List out(nb_channels);
    R_len_t tile_width = imgWidth / nb_channels;
    std::ifstream fi(fname.c_str(), std::ios::in | std::ios::binary | std::ios::ate);
    if (fi.is_open()) {
      try{
        fi.seekg(0, std::ios::end);
        unsigned int filesize = fi.tellg();
        if(verbose) {
          Rcout << fname << std::endl;
          Rcout << "Extracting " << nbytes << " Bytes BitMask image [30818] @offset:" << offset << std::endl;
        }
        if(offset > (filesize - nbytes)) {
          Rcpp::Rcerr << "hpp_rle_Decomp: @offset:" << offset << " points to outside of\n" << fname << std::endl;
          Rcpp::stop("hpp_rle_Decomp: RLE image offset is higher than file size");
        }
        fi.seekg(offset, std::ios::beg);
        std::vector<char> buf_image(nbytes);
        fi.read(buf_image.data(), nbytes);
        
        Rcpp::IntegerMatrix img(imgWidth,imgHeight);
        uint32_t L = imgWidth * imgHeight, runLength = 0;
        
        switch(removal) {
        case 1: { // clipped removal
          for(uint32_t k = 0; k < nbytes; k++) {
          int value = buf_image[k++];
          if(value > 1) value = -1;
          uint32_t off = runLength;
          runLength = off + (buf_image[k] & 0xff) + 1;
          if (runLength > L) {
            Rcpp::stop("hpp_rle_Decomp: Buffer overrun");
          }
          for(uint32_t j = off; j < runLength; j++) img[j] = value;
        }
          break;
        }
        case 2: { // height clipped removal
          for(uint32_t k = 0; k < nbytes; k++) {
          int value = buf_image[k++];
          if(value == 2) value = -1;
          uint32_t off = runLength;
          runLength = off + (buf_image[k] & 0xff) + 1;
          if (runLength > L) {
            Rcpp::stop("hpp_rle_Decomp: Buffer overrun");
          }
          for(uint32_t j = off; j < runLength; j++) img[j] = value;
        }
          break;
        }
        case 3: { // width clipped removal
          for(uint32_t k = 0; k < nbytes; k++) {
          int value = buf_image[k++];
          if(value == 3) value = -1;
          uint32_t off = runLength;
          runLength = off + (buf_image[k] & 0xff) + 1;
          if (runLength > L) {
            Rcpp::stop("hpp_rle_Decomp: Buffer overrun");
          }
          for(uint32_t j = off; j < runLength; j++) img[j] = value;
        }
          break;
        }
        case 4: { // only keep background
          for(uint32_t k = 0; k < nbytes; k++) {
          int value = buf_image[k++];
          value = (value == 0) ? 1:-1;
          uint32_t off = runLength;
          runLength = off + (buf_image[k] & 0xff) + 1;
          if (runLength > L) {
            Rcpp::stop("hpp_rle_Decomp: Buffer overrun");
          }
          for(uint32_t j = off; j < runLength; j++) img[j] = value;
        }
          break;
        }
        case 5: { // only keep non clipped foreground
          for(uint32_t k = 0; k < nbytes; k++) {
          int value = buf_image[k++];
          if(value != 1) value = -1;
          uint32_t off = runLength;
          runLength = off + (buf_image[k] & 0xff) + 1;
          if (runLength > L) {
            Rcpp::stop("hpp_rle_Decomp: Buffer overrun");
          }
          for(uint32_t j = off; j < runLength; j++) img[j] = value;
        }
          break;
        }
        default: { // no removal 
          for(uint32_t k = 0; k < nbytes; k++) {
          int value = buf_image[k++];
          uint32_t off = runLength;
          runLength = off + (buf_image[k] & 0xff) + 1;
          if (runLength > L) {
            Rcpp::stop("hpp_rle_Decomp: Buffer overrun");
          }
          for(uint32_t j = off; j < runLength; j++) img[j] = value;
        }
          break;
        }
        }
        fi.close();
        Rcpp::IntegerMatrix timg = Rcpp::transpose(img);
        for(R_len_t i = 0; i < nb_channels; i++) {
          out[i] = timg(Rcpp::_, Rcpp::Range(tile_width * i, tile_width * (i+1) - 1));
        }
        return out;
      }
      catch(std::exception &ex) {	
        fi.close();
        forward_exception_to_r(ex);
      }
      catch(...) { 
        Rcpp::stop("hpp_rle_Decomp: c++ exception (unknown reason)"); 
      }
    }
    else {
      Rcpp::stop("hpp_rle_Decomp: Unable to open file");
    }
  } else {
    Rcpp::stop("hpp_rle_Decomp: imgWidth, imgHeight and nb_channels should be >0");    
  }
  return R_NilValue;
}

//' @title GRAY Decompression
//' @name cpp_gray_Decomp
//' @description
//' Operates GrayScale decompression of compressed image stored in TIFF file.
//' @param fname string, path to file.
//' @param offset uint32_t, position of the beginning of compressed image.
//' @param nbytes uint32_t, number of bytes of compressed image.
//' @param imgWidth R_len_t, Width of the decompressed image. Default is 1.
//' @param imgHeight R_len_t, Height of the decompressed image. Default is 1.
//' @param nb_channels R_len_t, number of channels of the decompressed image. Default is 1.
//' @param verbose bool, whether to display information (use for debugging purpose). Default is false.
//' @details
//' BSD implementations of Bio-Formats readers and writers
//' %%
//' Copyright (C) 2005 - 2017 Open Microscopy Environment:
//'   - Board of Regents of the University of Wisconsin-Madison
//'   - Glencoe Software, Inc.
//'   - University of Dundee
//' %%
//' Redistribution and use in source and binary forms, with or without
//' modification, are permitted provided that the following conditions are met:
//' 
//' 1. Redistributions of source code must retain the above copyright notice,
//'    this list of conditions and the following disclaimer.
//' 2. Redistributions in binary form must reproduce the above copyright notice,
//'    this list of conditions and the following disclaimer in the documentation
//'    and/or other materials provided with the distribution.
//' 
//' THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
//' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
//' IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
//' ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
//' LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
//' CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
//' SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
//' INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
//' CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
//' ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
//' POSSIBILITY OF SUCH DAMAGE.
//' @source For image decompression, Lee Kamentsky's code porting from \url{https://github.com/ome/bioformats/blob/4146b9a1797501f0fec7d6cfe69124959bff96ee/components/formats-bsd/src/loci/formats/in/FlowSightReader.java}\cr
//' cited in \url{https://linkinghub.elsevier.com/retrieve/pii/S1046-2023(16)30291-2}
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::List hpp_gray_Decomp (const std::string fname, 
                            const uint32_t offset, 
                            const uint32_t nbytes,
                            const R_len_t imgWidth = 1, 
                            const R_len_t imgHeight = 1, 
                            const R_len_t nb_channels = 1,
                            const bool verbose = false) {
  if((nb_channels * imgWidth * imgHeight) != 0) {
    Rcpp::List out(nb_channels);
    R_len_t tile_width = imgWidth / nb_channels;
    std::ifstream fi(fname.c_str(), std::ios::in | std::ios::binary | std::ios::ate);
    if(fi.is_open()) {
      try {
        fi.seekg(0, std::ios::end);
        std::size_t filesize = fi.tellg();
        if(verbose) {
          Rcout << fname << std::endl;
          Rcout << "Extracting " << nbytes << " Bytes GreyScale image [30817] @offset:" << offset << std::endl;
        }
        if(offset > (filesize - nbytes)) {
          Rcpp::Rcerr << "hpp_gray_Decomp: @offset:" << offset << " points to outside of\n" << fname  << std::endl;
          Rcpp::stop("hpp_gray_Decomp: GrayScale image offset is higher than file size");
        }
        fi.seekg(offset, std::ios::beg);
        std::vector<char> buf_image(nbytes);
        fi.read(buf_image.data(), nbytes);
        
        Rcpp::IntegerVector lastRow(imgWidth + 1);
        Rcpp::IntegerMatrix img(imgHeight, imgWidth + 1);
        bool odd = false;
        
        uint32_t k = nbytes; // ensure that k will stay within [0, nbytes-1]
        for(R_len_t y = 0 ; y < imgHeight ; y++) {
          for(R_len_t x = 1 ; x <= imgWidth ; x++) {
            int value = 0;
            short shift = 0, nibble = -1;
            while((nibble & 0x8)) {
              nibble = odd ? buf_image[nbytes - k--] >> 4 : buf_image[nbytes - k] & 0xf;
              odd = !odd;
              value += (nibble & 0x7) << shift;
              shift += 3;
            }
            if(nibble & 0x4) value |= - (1 << shift);
            // img(y,x) = value;
            lastRow[x] += value;
            img(y,x) = img(y,x - 1) + lastRow[x];
          }
        }
        
        fi.close();
        for(R_len_t i = 0; i < nb_channels; i++) {
          out[i] = img(Rcpp::_, Rcpp::Range(1 + tile_width * i, tile_width * (i + 1)));
        }
        return out;
      }
      catch(std::exception &ex) {	
        fi.close();
        forward_exception_to_r(ex);
      }
      catch(...) { 
        Rcpp::stop("hpp_gray_Decomp: c++ exception (unknown reason)"); 
      }
    }
    else {
      Rcpp::stop("hpp_gray_Decomp: Unable to open file");
    }
  } else {
    Rcpp::stop("hpp_gray_Decomp: imgWidth, imgHeight and nb_channels should be >0");    
  }
  return R_NilValue;
}

//' @title IFC_object Decompression
//' @name cpp_decomp
//' @description
//' Operates decompression of compressed image stored in TIFF file.
//' @param fname string, path to file.
//' @param offset uint32_t, position of the beginning of compressed image.
//' @param nbytes uint32_t, number of bytes of compressed image.
//' @param imgWidth R_len_t, Width of the decompressed image. Default is 1.
//' @param imgHeight R_len_t, Height of the decompressed image. Default is 1.
//' @param nb_channels R_len_t, number of channels of the decompressed image. Default is 1.
//' @param removal uint8_t, object removal method. Only apply for 30818 compression. Default is 0.\cr
//' -1, for clipped removal: height OR width clipped pixels will be set to -1.\cr
//' -2, height clipped removal: height clipped pixels will be set to -1.\cr
//' -3, width clipped removal: width clipped pixels will be set to -1.\cr
//' -4, only keep background: background pixels will be set to 1 and all others to 0.\cr
//' -5, only keep foreground: foreground pixels will be set to 1 and all others to 0.
//' @param compression uint32_t, compression algorithm used. Default is 30818.
//' @param verbose bool, whether to display information (use for debugging purpose). Default is false.
//' @details
//' BSD implementations of Bio-Formats readers and writers
//' %%
//' Copyright (C) 2005 - 2017 Open Microscopy Environment:
//'   - Board of Regents of the University of Wisconsin-Madison
//'   - Glencoe Software, Inc.
//'   - University of Dundee
//' %%
//' Redistribution and use in source and binary forms, with or without
//' modification, are permitted provided that the following conditions are met:
//' 
//' 1. Redistributions of source code must retain the above copyright notice,
//'    this list of conditions and the following disclaimer.
//' 2. Redistributions in binary form must reproduce the above copyright notice,
//'    this list of conditions and the following disclaimer in the documentation
//'    and/or other materials provided with the distribution.
//' 
//' THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
//' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
//' IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
//' ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
//' LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
//' CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
//' SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
//' INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
//' CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
//' ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
//' POSSIBILITY OF SUCH DAMAGE.
//' @source For image decompression, Lee Kamentsky's code porting from \url{https://github.com/ome/bioformats/blob/4146b9a1797501f0fec7d6cfe69124959bff96ee/components/formats-bsd/src/loci/formats/in/FlowSightReader.java}\cr
//' cited in \url{https://linkinghub.elsevier.com/retrieve/pii/S1046-2023(16)30291-2}
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::List hpp_decomp (const std::string fname, 
                       const uint32_t offset, 
                       const uint32_t nbytes, 
                       const R_len_t imgWidth = 1, 
                       const R_len_t imgHeight = 1, 
                       const R_len_t nb_channels = 1,
                       const uint8_t removal = 0,
                       const uint32_t compression = 30818,
                       const bool verbose = false) {
  switch(compression) {
  case 30817: return hpp_gray_Decomp(fname, offset, nbytes, imgWidth, imgHeight, nb_channels, verbose);
  case 30818: return hpp_rle_Decomp(fname, offset, nbytes, imgWidth, imgHeight, nb_channels, removal, verbose);
  }
  Rcpp::Rcerr << "hpp_decomp: can't deal with compression format:" << compression << std::endl;
  Rcpp::stop("hpp_decomp: can't deal with compression format");   
  return R_NilValue;
}

#endif
