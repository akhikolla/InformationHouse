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

#ifndef IFC_TIFF_HPP
#define IFC_TIFF_HPP

#include <Rcpp.h>
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <iostream>
#include <fstream>
#include "utils.hpp"
using namespace Rcpp;

static int sizes[13] = {0,1,1,2,4,4,1,1,2,4,4,4,8};
static int multi[13] = {0,1,1,1,1,2,1,1,1,1,2,1,1};

// template to swap bytes for each types
template <typename T> T bytes_swap(T val) {
  T out;
  char *p_val = (char*)&val;
  char *p_out = (char*)&out;
  int size = sizeof(T);
  for(int i=0; i<size; i++) {
    p_out[size-1-i] = p_val[i];
  }
  return out;
}

//' @title TIFF Checker
//' @name cpp_checkTIFF
//' @description
//' Checks if file is a TIFF.
//' @details If file is a TIFF it returns endianness of file, 'big' or 'little.
//' Otherwise, it shows an error and returns an empty string.
//' @param fname string, path to file.
//' @source TIFF 6.0 specifications available at \url{https://www.adobe.io/open/standards/TIFF.html}
//' @keywords internal
////' @export
// [[Rcpp::export]]
std::string hpp_checkTIFF (const std::string fname) {
  std::ifstream fi(fname.c_str(), std::ios::in | std::ios::binary | std::ios::ate);
  std::string out = "";
  if (fi.is_open()) {
    try {
      fi.seekg(0, std::ios::end);
      std::size_t filesize = fi.tellg();
      if(filesize < 22) { // the 4 magic bytes + 1st Offset (4 bytes) + Count of Entry (2 bytes) + at least 1 IFD Entry (12 bytes) + last offset (4 bytes)
        Rcpp::stop("hpp_checkTIFF: File is too small");
      }
      uint16_t magic;
      fi.seekg(0, std::ios::beg);
      fi.read((char *)&magic,sizeof(magic));
      
      if(magic == 18761) out = "little"; // 49,49
      if(magic == 19789) out = "big"; // 4D,4D
      if(out == "") {
        Rcpp::stop("hpp_checkTIFF: File is not a XIF file: No magic bytes 0-1");
      }
      fi.read((char *)&magic,sizeof(magic));
      if(out == "big") magic = bytes_swap(magic);
      if(magic != 42) {
        Rcpp::stop("hpp_checkTIFF: File is not a XIF file: No magic bytes 2-3");
      }
      fi.close();
      return out;
    }
    catch(std::exception &ex) {	
      fi.close();
      forward_exception_to_r(ex);
    }
    catch(...) { 
      Rcpp::stop("hpp_checkTIFF: c++ exception (unknown reason)"); 
    }
  }
  else {
    Rcpp::stop("hpp_checkTIFF: Unable to open file");
  }
  return "";
}

//' @title IFC_offsets Computation without Id Determination
//' @name cpp_getoffsets_noid
//' @description
//' Returns offsets of the IFDs (Image Field Directory) within a TIFF file.
//' @param fname string, path to file.
//' @param obj_count R_len_t, numbers of objects present in the file. Default is 0.
//' If obj_count <= 0 then progress_bar is forced to false.
//' @param display_progress bool, whether to display a progress bar. Default is false.
//' @param verbose bool, whether to display information (use for debugging purpose). Default is false.
//' @source TIFF 6.0 specifications available at \url{https://www.adobe.io/open/standards/TIFF.html}
//' @return an integer vector with offsets of IFDs found.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::IntegerVector hpp_getoffsets_noid(const std::string fname, 
                                        const R_len_t obj_count = 0, 
                                        const bool display_progress = false,
                                        const bool verbose = false) {
  bool swap = false;
  std::string endianness = hpp_checkTIFF(fname);
  if(endianness == "big") swap = true;
  
  std::ifstream fi(fname.c_str(), std::ios::in | std::ios::binary | std::ios::ate);
  if (fi.is_open()) {
    try {
      fi.seekg(0, std::ios::end);
      std::size_t filesize = fi.tellg(); // offsets are uint32 so we can't go further
      uint32_t offset = 4; // offsets are uint32 
      std::size_t pos;
      Rcpp::IntegerVector out;
      bool show_pb = display_progress;
      if(obj_count <= 0) show_pb = false;
      char buf_entries [2];
      char buf_offset [4];
      uint16_t entries; // entries are uint16
      
      if(verbose) Rcout << "Extracting offsets from " << fname << std::endl;
      fi.seekg(offset, std::ios::beg);
      fi.read((char*)&buf_offset, sizeof(buf_offset));
      std::memcpy(&offset, buf_offset, sizeof(offset));
      
      Progress p(obj_count * 2 + 1, show_pb);
      if(swap) {
        offset = bytes_swap(offset);
        if(!offset) {
          Rcpp::stop("hpp_getoffsets_noid: No IFD offsets found");
        }
        while(offset) {
          p.increment();
          out.push_back(offset);
          fi.seekg(offset, std::ios::beg);
          fi.read((char*)buf_entries, sizeof(buf_entries));
          std::memcpy(&entries, buf_entries, sizeof(entries));
          entries = bytes_swap(entries);
          pos = 2 + offset + 12 * entries;
          if((pos > filesize) || (pos >= 0xffffffff)) { // can't read to more than 4,294,967,295
            Rcpp::stop("hpp_getoffsets_noid: Buffer overrun");
          }
          fi.seekg(pos, std::ios::beg);
          fi.read((char*)buf_offset, sizeof(buf_offset));
          std::memcpy(&offset, buf_offset, sizeof(offset));
          offset = bytes_swap(offset);
          if(Progress::check_abort()) {
            p.cleanup();
            Rcpp::stop("hpp_getoffsets_noid: Interrupted by user");
          }
        }
      } else {
        if(!offset) {
          Rcpp::stop("hpp_getoffsets_noid: No IFD offsets found");
        }
        while(offset) {
          p.increment();
          out.push_back(offset);
          fi.seekg(offset, std::ios::beg);
          fi.read((char*)buf_entries, sizeof(buf_entries));
          std::memcpy(&entries, buf_entries, sizeof(entries));
          pos = 2 + offset + 12 * entries;
          if((pos > filesize) || (pos >= 0xffffffff)) { // can't read to more than 4,294,967,295
            Rcpp::stop("hpp_getoffsets_noid: Buffer overrun");
          }
          fi.seekg(pos, std::ios::beg);
          fi.read((char*)buf_offset, sizeof(buf_offset));
          std::memcpy(&offset, buf_offset, sizeof(offset));
          if(Progress::check_abort()) {
            p.cleanup();
            Rcpp::stop("hpp_getoffsets_noid: Interrupted by user");
          }
        }
      }
      fi.close();
      return out;
    }
    catch(std::exception &ex) {	
      fi.close();
      forward_exception_to_r(ex);
    }
    catch(...) { 
      Rcpp::stop("hpp_getoffsets_noid: c++ exception (unknown reason)"); 
    }
  }
  else {
    Rcpp::stop("hpp_getoffsets_noid: Unable to open file");
  }
  return Rcpp::IntegerVector::create(0);
}

//' @title IFD Tags Extraction
//' @name cpp_getTAGS
//' @description
//' Returns TAGS contained within an IFD (Image Field Directory) entry.
//' @param fname string, path to file.
//' @param offset uint32_t, position of the IFD beginning.
//' @param verbose bool, whether to display information (use for debugging purpose). Default is 'false'.
//' @param trunc_bytes uint32_t maximal number of individual scalar to extract BYTE/ASCII/SBYTE/UNDIFINED for TAGS (1, 2, 6 or 7). Default is 12.\cr
//' However, if less is found, less is returned in map.
//' Note that, if 0 is provided, it will be automatically set to 1.
//' @param force_trunc whether to force truncation for all TAGS types. Default is FALSE.\cr
//' If 'true', 'trunc_bytes' will be used for TAGS (3, 4, 5, 8, 9, 10, 11 and 12) to extract desired number of individual scalar corresponding to each types.
//' @source TIFF 6.0 specifications available at \url{https://www.adobe.io/open/standards/TIFF.html}
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::List hpp_getTAGS (const std::string fname, 
                        const uint32_t offset, 
                        const bool verbose = false, 
                        const uint8_t trunc_bytes = 12, 
                        const bool force_trunc = false) {
  bool swap = false;
  std::string endianness = hpp_checkTIFF(fname); // it may be better to remove this check each time
  if(endianness == "big") swap = true; 
  
  std::ifstream fi(fname.c_str(), std::ios::in | std::ios::binary | std::ios::ate);
  if (fi.is_open()) {
    try {
      fi.seekg(0, std::ios::end);
      std::size_t filesize = fi.tellg();
      fi.seekg(0, std::ios::beg);
      if(offset > (filesize - 2)) {
        Rcpp::Rcerr << "hpp_getTAGS: offset@" << offset << " points to outside of " << fname << std::endl;
        Rcpp::stop("hpp_getTAGS: TAG offset is higher than file size");
      }
      if(verbose) Rcout << "Extracting TAGs from "<< fname << ", filesize:" << filesize << " @offset:" << offset << std::endl;
      
      uint16_t entries;
      uint16_t IFD_tag;
      uint16_t IFD_type;
      uint32_t IFD_count;
      uint32_t IFD_value;
      uint32_t IFD_bytes;
      uint32_t i;
      uint32_t tot_scalar;
      uint32_t ext_scalar;
      // ensure that, if not NULL, at least 1 scalar will be extracted
      uint32_t max_scalar = (trunc_bytes > 1) ? trunc_bytes:1;
      bool IFD_off;
      bool is_char;
      std::size_t pos;
      char buf_entries [2];
      char buf_dir_entry [12];
      
      fi.seekg(offset, std::ios::beg);
      fi.read((char*)&buf_entries, sizeof(buf_entries));
      std::memcpy(&entries, buf_entries, sizeof(entries));
      if(swap) entries = bytes_swap(entries);
      pos = offset + sizeof(buf_entries);
      
      RObject IMAGE_LENGTH = R_NilValue;          // 257
      RObject IMAGE_WIDTH = R_NilValue;           // 256
      // RObject TILE_LENGTH = R_NilValue;        // 323
      // RObject TILE_WIDTH = R_NilValue;         // 322
      RObject OBJECT_ID = R_NilValue;             // 33003
      // Rcpp::IntegerVector OBJECT_ID_2 = IR_NilValue; // 33024
      RObject COMPRESSION = R_NilValue;           // 259
      RObject TYPE = R_NilValue;                  // 33002
      RObject STRIP_OFFSETS = R_NilValue;         // 273
      RObject STRIP_BYTE_COUNTS = R_NilValue;     // 279
      // RObject TILE_OFFSETS = R_NilValue;       // 324
      // RObject TILE_BYTE_COUNTS = IR_NilValue;  // 325
      RObject BG_MEAN = R_NilValue;               // 33052
      RObject BG_STD = R_NilValue;                // 33053
      
      Rcpp::List INFOS(entries);
      CharacterVector NAMES(entries);
      if(verbose) Rcout << "Entries: " << entries << std::endl;
      if(swap) {
        for(uint16_t k = 0; k < entries; k++) {
          if((pos > (filesize - 12)) || (pos >= 0xffffffff)) {
            Rcpp::Rcerr << "hpp_getTAGS: in IFD: " << k << " @" << pos + 12 << " is outside of " << fname << std::endl;
            Rcpp::stop("hpp_getTAGS: IFD detected position is higher than file size");
          }
          fi.seekg(pos, std::ios::beg);
          fi.read((char*)&buf_dir_entry, sizeof(buf_dir_entry));
          std::memcpy(&IFD_tag, buf_dir_entry, sizeof(IFD_tag));
          IFD_tag = bytes_swap(IFD_tag);
          std::memcpy(&IFD_type, buf_dir_entry + 2, sizeof(IFD_type));
          IFD_type = bytes_swap(IFD_type);
          if((IFD_type > 12) || (IFD_type < 1)) {
            Rcpp::Rcerr << "hpp_getTAGS: in IFD: " << k << " IFD_type=" << IFD_type << " is not allowed" << std::endl;
            Rcpp::stop("hpp_getTAGS: Value not allowed for IFD type");
          }
          std::memcpy(&IFD_count, buf_dir_entry + 4, sizeof(IFD_count));
          IFD_count = bytes_swap(IFD_count);
          std::memcpy(&IFD_value, buf_dir_entry + 8, sizeof(IFD_value));
          IFD_value = bytes_swap(IFD_value);
          pos = fi.tellg();
          
          IFD_bytes = sizes[IFD_type] * multi[IFD_type] * IFD_count;
          tot_scalar = multi[IFD_type] * IFD_count;
          is_char = (IFD_type == 1) || (IFD_type == 2) || (IFD_type == 6) || (IFD_type == 7);
          
          // allow to extract less scalar than total number of scalars available
          // for char types or when force_trunc is choosen.
          // Except for tag 33052 (BG_MEAN) and 33053 (BG_STD)
          // that are retrieved completly (needed for img/msk extraction)
          ext_scalar = ((is_char || force_trunc) && (tot_scalar > max_scalar) && (IFD_tag != 33052) && (IFD_tag != 33053)) ? max_scalar:tot_scalar;
          if(IFD_bytes > 4) {
            if((IFD_value + IFD_bytes) > filesize) {
              Rcpp::Rcerr << "hpp_getTAGS: in IFD: " << k << " @" << IFD_value + IFD_bytes << " is outside of " << fname << std::endl;
              Rcpp::stop("hpp_getTAGS: IFD value points to outside of file");
            }
            fi.seekg(IFD_value, std::ios::beg);
            IFD_off = true;
          } else {
            fi.seekg(pos - 4, std::ios::beg);
            IFD_off = false;
          }
          if(verbose) Rcout << "Tag:" << IFD_tag << " Typ:" << IFD_type << " Count:" << IFD_count << " Value:" << IFD_value << " Bytes:" << IFD_bytes << " Off:" << IFD_off << std::endl;
          NAMES[k] = to_string(IFD_tag);
          
          if(is_char) {
            if(IFD_type == 2) {
              std::vector<char> buf_offset(ext_scalar * sizes[IFD_type]);
              fi.read(buf_offset.data(), ext_scalar * sizes[IFD_type]);
              std::string IFD_map(buf_offset.begin(), buf_offset.end());
              INFOS[k] = Rcpp::List::create(_["tag"] = IFD_tag,
                                            _["typ"] = IFD_type,
                                            _["siz"] = IFD_count,
                                            _["val"] = IFD_value,
                                            _["byt"] = IFD_bytes,
                                            _["len"] = tot_scalar,
                                            _["off"] = IFD_off,
                                            _["map"] = (IFD_count == 0) ? CharacterVector(0):IFD_map);
            } else {
              Rcpp::RawVector IFD_map(ext_scalar);
              switch(IFD_type) {
              case 1: {
                unsigned char ele;
                for(i = 0; i < ext_scalar; i++) {
                  fi.read((char *)&ele, sizes[IFD_type]);
                  IFD_map[i] = ele;
                }
              }
                break;
              case 6: {
                signed char ele;
                for(i = 0; i < ext_scalar; i++) {
                  fi.read((char *)&ele, sizes[IFD_type]);
                  IFD_map[i] = ele;
                }
              }
                break;
              case 7: {
                char ele;
                for(i = 0; i < ext_scalar; i++) {
                  fi.read((char *)&ele, sizes[IFD_type]);
                  IFD_map[i] = ele;
                }
              }
                break;
              }
              INFOS[k] = Rcpp::List::create(_["tag"] = IFD_tag,
                                            _["typ"] = IFD_type, 
                                            _["siz"] = IFD_count, 
                                            _["val"] = IFD_value, 
                                            _["byt"] = IFD_bytes, 
                                            _["len"] = tot_scalar, 
                                            _["off"] = IFD_off,
                                            _["map"] = IFD_map);
            }
          } else {
            Rcpp::NumericVector IFD_map(ext_scalar);
            switch(IFD_type) {
            case 3: {
              uint16_t ele;
              for(i = 0; i < ext_scalar; i++) {
                fi.read((char *)&ele, sizes[IFD_type]);
                IFD_map[i] = bytes_swap(ele);
              }
            }
              break;
            case 4: {
              uint32_t ele;
              for(i = 0; i < ext_scalar; i++) {
                fi.read((char *)&ele, sizes[IFD_type]);
                IFD_map[i] = bytes_swap(ele);
              }
            }
              break;
            case 5: {
              uint32_t ele;
              for(i = 0; i < ext_scalar; i++) {
                fi.read((char *)&ele, sizes[IFD_type]);
                IFD_map[i] = bytes_swap(ele);
              }
            }
              break;
            case 8: {
              int16_t ele;
              for(i = 0; i < ext_scalar; i++) {
                fi.read((char *)&ele, sizes[IFD_type]);
                IFD_map[i] = bytes_swap(ele);
              }
            }
              break;
            case 9: {
              int32_t ele;
              for(i = 0; i < ext_scalar; i++) {
                fi.read((char *)&ele, sizes[IFD_type]);
                IFD_map[i] = bytes_swap(ele);
              }
            }
              break;
            case 10: {
              int32_t ele;
              for(i = 0; i < ext_scalar; i++) {
                fi.read((char *)&ele, sizes[IFD_type]);
                IFD_map[i] = bytes_swap(ele);
              }
            }
              break;
            case 11: {
              float ele;
              for(i = 0; i < ext_scalar; i++) {
                fi.read((char *)&ele, sizes[IFD_type]);
                IFD_map[i] = bytes_swap(ele);
              }
            }
              break;
            case 12: {
              double ele;
              for(i = 0; i < ext_scalar; i++) {
                fi.read((char *)&ele, sizes[IFD_type]);
                IFD_map[i] = bytes_swap(ele);
              }
            }
              break;
            }
            
            switch(IFD_tag) {
            // TODO, ask AMNIS whether they use the TILE (OFFSETS/BYTE_COUNTS) or STRIP (OFFSETS/BYTE_COUNTS)
            // TODO, ask AMNIS whether they use the TILE (LENGTH/WIDTH) or IMAGE (LENGTH/WIDTH)
            case 256: {
              IMAGE_WIDTH = IFD_map[0];
              break;
            }
            case 257: {
              IMAGE_LENGTH = IFD_map[0];
              break;
            }
            case 259: {
              COMPRESSION = IFD_map[0];
              break;
            }
            case 273: {
              STRIP_OFFSETS = IFD_map[0];
              break;
            }
            case 279: {
              STRIP_BYTE_COUNTS = IFD_map[0];
              break;
            }
            case 33002: {
              TYPE = IFD_map[0];
              break;
            }
            case 33003: {
              OBJECT_ID = IFD_map[0];
              break;
            }
            case 33052: {
              BG_MEAN = IFD_map;
              break;
            }
            case 33053: {
              BG_STD = IFD_map;
              break;
            }
            }
            INFOS[k] = Rcpp::List::create(_["tag"] = IFD_tag,
                                          _["typ"] = IFD_type, 
                                          _["siz"] = IFD_count, 
                                          _["val"] = IFD_value, 
                                          _["byt"] = IFD_bytes, 
                                          _["len"] = tot_scalar, 
                                          _["off"] = IFD_off,
                                          _["map"] = IFD_map);
          }
        }
      } else {
        for(uint16_t k = 0; k < entries; k++) {
          if((pos > (filesize - 12)) || (pos >= 0xffffffff)) {
            Rcpp::Rcerr << "hpp_getTAGS: in IFD: " << k << " @" << pos + 12 << " is outside of " << fname << std::endl;
            Rcpp::stop("hpp_getTAGS: IFD detected position is higher than file size");
          }
          fi.seekg(pos, std::ios::beg);
          fi.read((char*)&buf_dir_entry, sizeof(buf_dir_entry));
          std::memcpy(&IFD_tag, buf_dir_entry, sizeof(IFD_tag));
          std::memcpy(&IFD_type, buf_dir_entry + 2, sizeof(IFD_type));
          if((IFD_type > 12) || (IFD_type < 1)) {
            Rcpp::Rcerr << "hpp_getTAGS: in IFD: " << k << " IFD_type=" << IFD_type << " is not allowed" << std::endl;
            Rcpp::stop("hpp_getTAGS: Value not allowed for IFD type");
          }
          std::memcpy(&IFD_count, buf_dir_entry + 4, sizeof(IFD_count));
          std::memcpy(&IFD_value, buf_dir_entry + 8, sizeof(IFD_value));
          pos = fi.tellg();
          
          IFD_bytes = sizes[IFD_type] * multi[IFD_type] * IFD_count;
          tot_scalar = multi[IFD_type] * IFD_count;
          is_char = (IFD_type == 1) || (IFD_type == 2) || (IFD_type == 6) || (IFD_type == 7);
          
          // allow to extract less scalar than total number of scalars available
          // for char types or when force_trunc is choosen.
          // Except for tag 33052 (BG_MEAN) and 33053 (BG_STD)
          // that are retrieved completly (needed for img/msk extraction)
          ext_scalar = ((is_char || force_trunc) && (tot_scalar > max_scalar) && (IFD_tag != 33052) && (IFD_tag != 33053)) ? max_scalar:tot_scalar;
          if(IFD_bytes > 4) {
            if((IFD_value + IFD_bytes) > filesize) {
              Rcpp::Rcerr << "hpp_getTAGS: in IFD: " << k << " @" << IFD_value + IFD_bytes << " is outside of " << fname << std::endl;
              Rcpp::stop("hpp_getTAGS: IFD value points to outside of file");
            }
            fi.seekg(IFD_value, std::ios::beg);
            IFD_off = true;
          } else {
            fi.seekg(pos - 4, std::ios::beg);
            IFD_off = false;
          }
          if(verbose) Rcout << "Tag:" << IFD_tag << " Typ:" << IFD_type << " Count:" << IFD_count << " Value:" << IFD_value << " Bytes:" << IFD_bytes << " Off:" << IFD_off << std::endl;
          NAMES[k] = to_string(IFD_tag);
          
          if(is_char) {
            if(IFD_type == 2) {
              std::vector<char> buf_offset(ext_scalar * sizes[IFD_type]);
              fi.read(buf_offset.data(), ext_scalar * sizes[IFD_type]);
              std::string IFD_map(buf_offset.begin(), buf_offset.end());
              INFOS[k] = Rcpp::List::create(_["tag"] = IFD_tag,
                                            _["typ"] = IFD_type,
                                            _["siz"] = IFD_count,
                                            _["val"] = IFD_value,
                                            _["byt"] = IFD_bytes,
                                            _["len"] = tot_scalar,
                                            _["off"] = IFD_off,
                                            _["map"] = (IFD_count == 0) ? CharacterVector(0):IFD_map);
            } else {
              Rcpp::RawVector IFD_map(ext_scalar);
              switch(IFD_type) {
              case 1: {
                unsigned char ele;
                for(i = 0; i < ext_scalar; i++) {
                  fi.read((char *)&ele, sizes[IFD_type]);
                  IFD_map[i] = ele;
                }
              }
                break;
              case 6: {
                signed char ele;
                for(i = 0; i < ext_scalar; i++) {
                  fi.read((char *)&ele, sizes[IFD_type]);
                  IFD_map[i] = ele;
                }
              }
                break;
              case 7: {
                char ele;
                for(i = 0; i < ext_scalar; i++) {
                  fi.read((char *)&ele, sizes[IFD_type]);
                  IFD_map[i] = ele;
                }
              }
                break;
              }
              INFOS[k] = Rcpp::List::create(_["tag"] = IFD_tag,
                                            _["typ"] = IFD_type, 
                                            _["siz"] = IFD_count, 
                                            _["val"] = IFD_value, 
                                            _["byt"] = IFD_bytes, 
                                            _["len"] = tot_scalar, 
                                            _["off"] = IFD_off,
                                            _["map"] = IFD_map);
            }
          } else {
            Rcpp::NumericVector IFD_map(ext_scalar);
            switch(IFD_type) {
            case 3: {
              uint16_t ele;
              for(i = 0; i < ext_scalar; i++) {
                fi.read((char *)&ele, sizes[IFD_type]);
                IFD_map[i] = ele;
              }
            }
              break;
            case 4: {
              uint32_t ele;
              for(i = 0; i < ext_scalar; i++) {
                fi.read((char *)&ele, sizes[IFD_type]);
                IFD_map[i] = ele;
              }
            }
              break;
            case 5: {
              uint32_t ele;
              for(i = 0; i < ext_scalar; i++) {
                fi.read((char *)&ele, sizes[IFD_type]);
                IFD_map[i] = ele;
              }
            }
              break;
            case 8: {
              int16_t ele;
              for(i = 0; i < ext_scalar; i++) {
                fi.read((char *)&ele, sizes[IFD_type]);
                IFD_map[i] = ele;
              }
            }
              break;
            case 9: {
              int32_t ele;
              for(i = 0; i < ext_scalar; i++) {
                fi.read((char *)&ele, sizes[IFD_type]);
                IFD_map[i] = ele;
              }
            }
              break;
            case 10: {
              int32_t ele;
              for(i = 0; i < ext_scalar; i++) {
                fi.read((char *)&ele, sizes[IFD_type]);
                IFD_map[i] = ele;
              }
            }
              break;
            case 11: {
              float ele;
              for(i = 0; i < ext_scalar; i++) {
                fi.read((char *)&ele, sizes[IFD_type]);
                IFD_map[i] = ele;
              }
            }
              break;
            case 12: {
              double ele;
              for(i = 0; i < ext_scalar; i++) {
                fi.read((char *)&ele, sizes[IFD_type]);
                IFD_map[i] = ele;
              }
            }
              break;
            }
            
            switch(IFD_tag) {
            // TODO, ask AMNIS whether they use the TILE (OFFSETS/BYTE_COUNTS) or STRIP (OFFSETS/BYTE_COUNTS)
            // TODO, ask AMNIS whether they use the TILE (LENGTH/WIDTH) or IMAGE (LENGTH/WIDTH)
            case 256: {
              IMAGE_WIDTH = IFD_map[0];
              break;
            }
            case 257: {
              IMAGE_LENGTH = IFD_map[0];
              break;
            }
            case 259: {
              COMPRESSION = IFD_map[0];
              break;
            }
            case 273: {
              STRIP_OFFSETS = IFD_map[0];
              break;
            }
            case 279: {
              STRIP_BYTE_COUNTS = IFD_map[0];
              break;
            }
            case 33002: {
              TYPE = IFD_map[0];
              break;
            }
            case 33003: {
              OBJECT_ID = IFD_map[0];
              break;
            }
            case 33052: {
              BG_MEAN = IFD_map;
              break;
            }
            case 33053: {
              BG_STD = IFD_map;
              break;
            }
            }
            INFOS[k] = Rcpp::List::create(_["tag"] = IFD_tag,
                                          _["typ"] = IFD_type, 
                                          _["siz"] = IFD_count, 
                                          _["val"] = IFD_value, 
                                          _["byt"] = IFD_bytes, 
                                          _["len"] = tot_scalar, 
                                          _["off"] = IFD_off,
                                          _["map"] = IFD_map);
          }
        }
      }
      
      char buf_next [4];
      unsigned int next;
      fi.seekg(pos, std::ios::beg);
      fi.read((char*)&buf_next, sizeof(buf_next));
      std::memcpy(&next, buf_next, sizeof(next));
      if(swap) next = bytes_swap(next); 
      INFOS.names() = NAMES;
      Rcpp::List out = Rcpp::List::create(_["tags"] = INFOS,
                                          _["infos"] = Rcpp::List::create(
                                            _["IMAGE_LENGTH"] = IMAGE_LENGTH,
                                            _["IMAGE_WIDTH"] = IMAGE_WIDTH,
                                            // _["TILE_LENGTH"] = TILE_LENGTH,
                                            // _["TILE_WIDTH"] = TILE_WIDTH,
                                            _["OBJECT_ID"] = OBJECT_ID,
                                            // _["OBJECT_ID_2"] = OBJECT_ID_2,
                                            _["COMPRESSION"] = COMPRESSION,
                                            _["TYPE"] = TYPE,
                                            _["STRIP_OFFSETS"] = STRIP_OFFSETS,
                                            _["STRIP_BYTE_COUNTS"] = STRIP_BYTE_COUNTS,
                                            // _["TILE_OFFSETS"] = TILE_OFFSETS,
                                            // _["TILE_BYTE_COUNTS"] = TILE_BYTE_COUNTS,
                                            _["BG_MEAN"] = BG_MEAN,
                                            _["BG_STD"] = BG_STD),
                                            _["curr_IFD_offset"] = offset,
                                            _["next_IFD_offset"] = next);
      out.attr("class") = "IFC_ifd";
      fi.close();
      return out;
    }
    catch(std::exception &ex) {
      fi.close();
      forward_exception_to_r(ex);
    }
    catch(...) { 
      Rcpp::stop("hpp_getTAGS: c++ exception (unknown reason)"); 
    }
  }
  else {
    Rcpp::stop("hpp_getTAGS: Unable to open file");
  }
  return Rcpp::List::create(_["tags"] = NA_REAL,
                            _["infos"] = Rcpp::List::create(
                              _["IMAGE_LENGTH"] = NA_REAL,
                              _["IMAGE_WIDTH"] = NA_REAL,
                              // _["TILE_LENGTH"] = NA_REAL,
                              // _["TILE_WIDTH"] = NA_REAL,
                              _["OBJECT_ID"] = NA_REAL,
                              // _["OBJECT_ID_2"] = NA_REAL,
                              _["COMPRESSION"] = NA_REAL,
                              _["TYPE"] = NA_REAL,
                              _["STRIP_OFFSETS"] = NA_REAL,
                              _["STRIP_BYTE_COUNTS"] = NA_REAL,
                              // _["TILE_OFFSETS"] = NA_REAL,
                              // _["TILE_BYTE_COUNTS"] = NA_REAL,
                              _["BG_MEAN"] = NA_REAL,
                              _["BG_STD"] = NA_REAL),
                              _["curr_IFD_offset"] = NA_REAL,
                              _["next_IFD_offset"] = NA_REAL);
}

//' @title IFC_offsets Computation with Object Identification
//' @name cpp_getoffsets_wid
//' @description
//' Returns offsets of the IFD (Image Field Directory) within a TIFF file.
//' @param fname string, path to file.
//' @param obj_count R_len_t, numbers of objects present in the file. Default is 0.
//' If obj_count <= 0 then progress_bar is forced to false.
//' @param display_progress bool, whether to display a progress bar. Default is false.
//' @param verbose bool, whether to display information (use for debugging purpose). Default is false.
//' @source TIFF 6.0 specifications available at \url{https://www.adobe.io/open/standards/TIFF.html}
//' @return a list of integer vectors with OBJECT_ID, TYPE and OFFSET of IFDs found.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::List hpp_getoffsets_wid(const std::string fname, 
                              const R_len_t obj_count = 0, 
                              const bool display_progress = false, 
                              const bool verbose = false) {
  bool swap = false;
  std::string endianness = hpp_checkTIFF(fname);
  if(endianness == "big") swap = true;
  
  std::ifstream fi(fname.c_str(), std::ios::in | std::ios::binary | std::ios::ate);
  if (fi.is_open()) {
    try{
      fi.seekg(0, std::ios::end);
      bool show_pb = display_progress;
      if(obj_count == 0) show_pb = false;
      char buf_offset [4];
      uint32_t offset = 4;
      
      if(verbose) Rcout << "Extracting offsets from " << fname << std::endl;
      Progress p(obj_count * 2 + 1, show_pb);
      fi.seekg(offset, std::ios::beg);
      fi.read((char*)&buf_offset, sizeof(buf_offset));
      std::memcpy(&offset, buf_offset, sizeof(offset));
      if(swap) offset = bytes_swap(offset);
      if(!offset) {
        Rcpp::stop("hpp_getoffsets_wid: No IFD offsets found");
      }
      
      Rcpp::IntegerVector out_obj, out_typ, out_off;
      // Rcpp::IntegerVector obj = Rcpp::IntegerVector::create(NA_INTEGER);
      // Rcpp::IntegerVector typ = Rcpp::IntegerVector::create(NA_INTEGER);
      // uint32_t obj, typ;
      
      // uint32_t trunc_bytes = 8;
      // bool force_trunc = true;
      while(offset){
        p.increment();
        Rcpp::List IFD = hpp_getTAGS(fname, offset, verbose, 4, true);
        offset = IFD["next_IFD_offset"];
        Rcpp::List infos = IFD["infos"];
        
        if(iNotisNULL(infos["OBJECT_ID"])) {
          out_obj.push_back(infos["OBJECT_ID"]);
          // obj[0] = clone(std::as<int32_t>(infos["OBJECT_ID"])); // min 0?, max ?
        } else {
          out_obj.push_back(NA_INTEGER);
        }
        if(iNotisNULL(infos["TYPE"])) {
          out_typ.push_back(infos["TYPE"]);
          // typ[0] = clone(std::as<int32_t>(infos["TYPE"])); // [1,3]
        } else {
          out_typ.push_back(NA_INTEGER);
        }
        
        // out_obj.push_back(obj[0]);
        // out_typ.push_back(typ[0]);
        out_off.push_back(IFD["curr_IFD_offset"]);
        
        if(Progress::check_abort()) {
          p.cleanup();
          Rcpp::stop("hpp_getoffsets_wid: Interrupted by user");
        }
      }
      Rcpp::List out = Rcpp::List::create(_["OBJECT_ID"] = out_obj,
                                          _["TYPE"] = out_typ,
                                          _["OFFSET"] = out_off);
      fi.close();
      return out;
    }
    catch(std::exception &ex) {	
      fi.close();
      forward_exception_to_r(ex);
    }
    catch(...) { 
      Rcpp::stop("hpp_getoffsets_wid: c++ exception (unknown reason)"); 
    }
  }
  else {
    Rcpp::stop("hpp_getoffsets_wid: Unable to open file");
  }
  return Rcpp::List::create(Rcpp::List::create(_["OBJECT_ID"] = NA_INTEGER,
                                               _["TYPE"] = NA_INTEGER,
                                               _["OFFSET"] = NA_INTEGER));
}

//' @title Checksum for RIF/CIF
//' @name cpp_checksum
//' @description
//' Computes sum of img IFDs (Image Field Directory) offsets of objects 0, 1, 2, 3 and 4.
//' @param fname string, path to file.
//' @source TIFF 6.0 specifications available at \url{https://www.adobe.io/open/standards/TIFF.html}
//' @return an integer vector with offsets of IFDs found.
//' @keywords internal
////' @export
// [[Rcpp::export]]
std::size_t hpp_checksum(const std::string fname) {
  bool swap = false;
  std::string endianness = hpp_checkTIFF(fname);
  if(endianness == "big") swap = true;
  
  std::ifstream fi(fname.c_str(), std::ios::in | std::ios::binary | std::ios::ate);
  if (fi.is_open()) {
    try {
      fi.seekg(0, std::ios::end);
      char buf_offset [4];
      uint32_t offset = 4;
      Rcpp::IntegerVector obj = Rcpp::IntegerVector::create(0,1,2,3,4);
      Rcpp::IntegerVector found = Rcpp::IntegerVector::create(6); // ensure to test at least one value in found == id
      uint8_t count = 0;
      std::size_t out = 0;
      
      fi.seekg(offset, std::ios::beg);
      fi.read((char*)&buf_offset, sizeof(buf_offset));
      std::memcpy(&offset, buf_offset, sizeof(offset));
      if(swap) offset = bytes_swap(offset);
      if(!offset) {
        Rcpp::stop("hpp_checksum: No IFD offsets found");
      }
      bool warn = true;
      while(offset && (count < 5)){
        Rcpp::List IFD = hpp_getTAGS(fname, offset, false, 8, true);
        Rcpp::List infos = IFD["infos"];
        offset = as<uint32_t>(IFD["next_IFD_offset"]);
        int32_t typ = 0;
        // infos["TYPE"] should never be NULL
        if(iNotisNULL(infos["TYPE"])) typ = as<int32_t>(infos["TYPE"]); // [1,3]
        // infos["OBJECT_ID"] should be NULL when typ is 1 or 3 but never when typ is 2
        if(iNotisNULL(infos["OBJECT_ID"]) && (typ==2)) { // ensure it is image
          int32_t id = as<int32_t>(infos["OBJECT_ID"]);
          if(is_true(any(obj == id)) && !is_true(any(found == id))) {
            found.push_back(id);
            out += as<uint32_t>(IFD["curr_IFD_offset"]);
            if(warn) if(id != obj[count]) { // ensure it is stored in ascending order
              Rcpp::warning("hpp_checksum: raw object are not stored in expected order");
              warn = false;
            }
            count += 1; 
          } else { // ensure it is stored in ascending order
            if(warn) {
              Rcpp::warning("hpp_checksum: raw object are not stored in expected order");
              warn = false;
            }
          }
        }
      }
      fi.close();
      return out;
    }
    catch(std::exception &ex) {	
      fi.close();
      forward_exception_to_r(ex);
    }
    catch(...) { 
      Rcpp::stop("hpp_checksum: c++ exception (unknown reason)"); 
    }
  }
  else {
    Rcpp::stop("hpp_checksum: Unable to open file");
  }
  return 0;
}

#endif
