################################################################################
# This file is released under the GNU General Public License, Version 3, GPL-3 #
# Copyright (C) 2020 Yohann Demont                                             #
#                                                                              #
# It is part of IFC package, please cite:                                      #
# -IFC: An R Package for Imaging Flow Cytometry                                #
# -YEAR: 2020                                                                  #
# -COPYRIGHT HOLDERS: Yohann Demont, Gautier Stoll, Guido Kroemer,             #
#                     Jean-Pierre Marolleau, Loïc Garçon,                      #
#                     INSERM, UPD, CHU Amiens                                  #
#                                                                              #
# DISCLAIMER:                                                                  #
# -You are using this package on your own risk!                                #
# -We do not guarantee privacy nor confidentiality.                            #
# -This program is distributed in the hope that it will be useful, but WITHOUT #
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or        #
# FITNESS FOR A PARTICULAR PURPOSE. In no event shall the copyright holders or #
# contributors be liable for any direct, indirect, incidental, special,        #
# exemplary, or consequential damages (including, but not limited to,          #
# procurement of substitute goods or services; loss of use, data, or profits;  #
# or business interruption) however caused and on any theory of liability,     #
# whether in contract, strict liability, or tort (including negligence or      #
# otherwise) arising in any way out of the use of this software, even if       #
# advised of the possibility of such damage.                                   #
#                                                                              #
# You should have received a copy of the GNU General Public License            #
# along with IFC. If not, see <http://www.gnu.org/licenses/>.                  #
################################################################################

#' @title Image Field Directory Builder
#' @description Builds Image Field Directory (IFD)
#' @param val the value of the IFD
#' @param typ desired IFD type
#' @param tag the desired IFD 'tag'
#' @param endianness the desired endian-ness ("big" or "little"). Default is .Platform$endian.\cr
#' Endianness describes the bytes order of data stored within the files. This parameter may not be modified.
#' @details if 'val' if of type "character", 'tag' is automatically set to 2.\cr
#' if 'val' is of length 0 NULL is returned.
#' @return NULL or a list of 2 members:\cr
#' -min_content: the minimal IFD content,\cr
#' -add_content: the additional IFD content if 'val' converted to raw does not fit in 4 bytes.
#' @keywords internal
buildIFD <- function(val, typ, tag, endianness = .Platform$endian) {
  sizes = c(1,1,2,4,4,1,1,2,4,4,4,8)
  multi = c(1,1,1,1,2,1,1,1,1,2,1,1)
  switch(typeof(val),
         "character" = { 
           typ <- 2
           val = strsplit(x = unname(val), split = character())
           if(length(val) == 1) val = val[[1]]
         })
  val_raw = lapply(unname(val), FUN = function(x) {
    switch(typ,
           { x # 1 BYTE
           },
           { charToRaw(as.character(x)) # 2 ASCII
           },
           { packBits(intToBits(x),type="raw")[1:2] # 3 SHORT 2 bytes, what happen when endianness is swapped ?
           },
           { packBits(intToBits(x),type="raw") # 4 LONG, 4 bytes
           },
           { packBits(intToBits(x),type="raw") # 5 RATIONAL = 2 LONG
           },
           { x # 6 SBYTE
           },
           { x # 7 UNDEFINED, 1 Byte
           },
           { packBits(intToBits(x),type="raw")[1:2] # 8 SSHORT, 2 bytes, what happen when endianness is swapped ?
           },
           { packBits(intToBits(x),type="raw") # 9 SLONG, 4 bytes
           },
           { packBits(intToBits(x),type="raw") # 10 SRATIONAL, 2 SLONG
           },
           { writeBin(x, raw(), size = 4) # 11 FLOAT, 4 bytes
           },
           { writeBin(x, raw(), size = 8) # 12 DOUBLE, 8 bytes
           })
  })
  bytes = length(unlist(val_raw))
  count = bytes / (sizes[typ] * multi[typ])
  ifd = list(packBits(intToBits(tag),type="raw")[1:2], #tag
             packBits(intToBits(typ),type="raw")[1:2], #typ
             packBits(intToBits(count),type="raw")) #count
  if(endianness != .Platform$endian) {
    ifd = lapply(ifd, rev)
    val_raw = lapply(val_raw, rev)
  }
  if(bytes > 4) {
    ifd = c(ifd, packBits(intToBits(0),type="raw")) #val/offsets
    add = val_raw
  } else {
    ifd = c(ifd, sapply(1:4, FUN=function(i) { ifelse(i <= bytes, unlist(val_raw)[i], as.raw(0x00)) }))
    add = raw()
  }
  ifd = list(list(min_content = as.raw(unlist(ifd)), add_content = unlist(add)))
  names(ifd) = tag
  return(ifd)
}
