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

#ifndef IFC_SCAN_HPP
#define IFC_SCAN_HPP

#include <Rcpp.h>
#include <iostream>
#include <fstream>
using namespace Rcpp;

//' @title File Scanner
//' @name cpp_scanFirst
//' @description
//' Scans file for 1st occurence of a target string.
//' If found, it returns the position in bytes of the target.
//' Otherwise, it returns 0.
//' @param fname string, path to file.
//' @param target string, exact string to be searched for. At least 1 character and should not exceed 1024 characters.
//' @param start size_t, position where to begin search.
//' It can't be superior or equal than file size or end (when end is different from 0 and inferior than file size).
//' @param end size_t, position where to stop searching. Default is 0.
//' Search will end up at this position unless it is higher than file size.
//' In such case, search will end up when file end will be reached.
//' @param buf_size uint8_t, size of buffer used to search for target (in kilo-Bytes, will be forced to be between 2 and 1024). Default is 64.
//' @return size_t index of first target character found within target plus 1 or 0 if not found.
//' @keywords internal
////' @export
// [[Rcpp::export]]
std::size_t hpp_scanFirst(const std::string fname, 
                          const std::string target, 
                          const std::size_t start = 0, 
                          const std::size_t end = 0, 
                          const uint8_t buf_size = 64) {
  uint16_t L = target.length();
  if(L < 1) {
    Rcpp::stop("cpp_scanFirst: target should be at least 1 character");
  }
  if(L > 1024) {
    Rcpp::Rcerr <<  "cpp_scanFirst: target should not exceed 1024 characters but is of length: [" << L << "], for file:" << std::endl << fname << std::endl;
    Rcpp::stop("cpp_scanFirst: target should not exceed 1024 characters");
  }
  
  std::ifstream fi(fname.c_str(), std::ios::in | std::ios::binary | std::ios::ate);
  if (fi.is_open()) {
    try {
      fi.seekg(0, std::ios::end);
      std::size_t filesize = fi.tellg(), end_at = 0, n = 0, pos = 0;
      bool keep_searching = true;
      
      // determines where to stop
      end_at = filesize;
      if(end != 0) if(end < filesize) end_at = end;
      if(start > (end_at - L)) return 0;
      
      // clip buffer to [2,255]
      uint8_t buf_s = buf_size;
      if(buf_s < 2) buf_s = 2;
      uint32_t s = 1024 * buf_s; // allow reading chunk of size buf_s * kilo-Bytes
      
      fi.seekg(start, std::ios::beg);
      while(keep_searching) {
        pos = fi.tellg();
        if(pos >= end_at) break; // break loop when end of the file is reached (or end when provided);
        if(pos > (start + L)) pos -= L; // make a frame shift for every chunk except 1st one to be sure to include target when it is positioned between 2 chunks;
        fi.seekg(pos, std::ios::beg);
        if((end_at - pos) < s) s = end_at - pos; // reduce buffer size for last chunk to avoid reading outside of the file;
        std::vector<char> buffer(s);
        fi.read(buffer.data(), s);
        std::string str(buffer.begin(), buffer.end());
        n = str.find(target);
        if(n != std::string::npos) keep_searching = false;
      }
      fi.close();
      if(!keep_searching) return pos + n + 1;
    }
    catch(std::exception &ex) { 
      fi.close();
      forward_exception_to_r(ex);
    }
    catch(...) {
      Rcpp::stop("cpp_scanFirst: c++ exception (unknown reason)");
    }
  }
  else {
    Rcpp::stop("cpp_scanFirst: Unable to open file");
  }
  return 0;
}

#endif
