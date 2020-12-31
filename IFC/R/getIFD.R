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

#' @title RIF/CIF Image Field Directories Extraction
#' @description
#' Extracts IFDs (Image File Directory) in RIF or CIF files.\cr
#' IFDs contain information about images or masks of objects stored within XIF files.\cr
#' The first IFD is special in that it does not contain image of mask information but general information about the file.\cr
#' Users are highly encouraged to read TIFF specifications to have a better understanding about what IFDs are.
#' @param fileName path to file.
#' @param offsets either "all", "first" or an object of class `IFC_offset`. Default is "first".
#' @param trunc_bytes a positive integer maximal number of individual scalar to extract BYTE/ASCII/SBYTE/UNDIFINED for TAGS (1, 2, 6 or 7). Default is 12.\cr
#' However, if less is found, less is returned in map. 
#' Note that, if 0 is provided, it will be automatically set to 1.
#' @param force_trunc whether to force truncation for all TAGS types. Default is FALSE.\cr
#' If TRUE, 'trunc_bytes' will be used for TAGS (3, 4, 5, 8, 9, 10, 11 and 12) to extract desired number of individual scalar corresponding to each types.
#' @param verbose whether to display information (use for debugging purpose). Default is FALSE.
#' @param verbosity quantity of information displayed when verbose is TRUE; 1: normal, 2: rich. Default is 1.
#' @param display_progress whether to display a progress bar. Default is FALSE.
#' @param bypass whether to bypass checks on 'trunc_bytes', 'force_trunc', 'verbose', 'verbosity' and 'display_progress'. Default is FALSE.
#' @param ... other arguments to be passed.
#' @source TIFF 6.0 specifications available at \url{https://www.adobe.io/open/standards/TIFF.html}
#' @details Function will return IFDs (image, mask or first) from the file using provided offsets argument.\cr
#' IFDs contain several tags that can be viewed as descriptive meta-information of raw data stored within RIF or CIF file. For more details see TIFF specifications.\cr
#' If 'offsets' == "first" only first IFD will be returned.\cr
#' If 'offsets' == "all" all images and masks IFDs will be returned but not "first" one.
#' Be aware that errors may occur if offsets are not extracted with \code{\link{getOffsets}} or \code{\link{subsetOffsets}}.
#' @examples
#' if(requireNamespace("IFCdata", quietly = TRUE)) {
#'   ## use a cif file
#'   file_cif <- system.file("extdata", "example.cif", package = "IFCdata")
#'   ## read 1st IFD
#'   IFD_first <- getIFD(fileName = file_cif, offsets = "first")
#'   ## show information contained in 1st IFD
#'   print(sapply(IFD_first[[1]]$tags, FUN=function(x) x)) 
#' } else {
#'   message(sprintf('Please run `install.packages("IFCdata", repos = "%s", type = "source")` %s',
#'                   'https://gitdemont.github.io/IFCdata/',
#'                   'to install extra files required to run this example.'))
#' }
#' @return A list of named lists, each containing:\cr
#' -tags, a named list whose names are tags found, where each tag is a list of tag, typ, siz, val, byt, len, off, map information.\cr
#' -infos, a named list containing essential information about IFDs,
#' IMAGE_LENGTH, IMAGE_WIDTH, OBJECT_ID, COMPRESSION,
#' TYPE, STRIP_OFFSETS, STRIP_BYTE_COUNTS,
#' BG_MEAN, BG_STD\cr
#' -curr_IFD_offset, the position of current IFD offset\cr
#' -next_IFD_offset, the position of next IFD offset
#' @export
getIFD <- function(fileName, offsets = "first", trunc_bytes = 12, force_trunc = FALSE,
                   verbose = FALSE, verbosity = 1, display_progress = FALSE, bypass = FALSE, ...) {
  dots = list(...)
  # various check
  assert(bypass, len = 1, alw = c(TRUE, FALSE))
  endianness = cpp_checkTIFF(fileName) # used to determine endianness and check that file exists and is of XIF content.
  title_progress = basename(fileName)
  if(!bypass) {  # bypass checking
    trunc_bytes = as.integer(trunc_bytes); assert(trunc_bytes, len = 1, typ = "integer")
    force_trunc = as.logical(force_trunc); assert(force_trunc, len = 1, alw = c(TRUE, FALSE))
    verbose = as.logical(verbose); assert(verbose, len = 1, alw = c(TRUE, FALSE))
    verbosity = as.integer(verbosity); assert(verbosity, len = 1, alw = c(1, 2))
    display_progress = as.logical(display_progress); assert(display_progress, len = 1, alw = c(TRUE, FALSE))
  }
  
  # open fileName to extract 1st offset
  toread = file(fileName,"rb")
  on.exit(close(toread), add = TRUE)
  first_offset = 4
  seek(toread, first_offset)
  first_offset = readBin(toread, "integer", n = 1, endian = endianness)
  # extract 1st IFD
  first_IFD = cpp_getTAGS(fname = fileName, offset = first_offset, trunc_bytes = trunc_bytes, force_trunc = force_trunc, verbose = ifelse(verbose & (verbosity==2), TRUE, FALSE))
  
  # extract important information from 1st IFD
  obj_number = first_IFD$tags$`33018`$map

  # checks offsets 
  if(length(offsets) == 1) {
    if(offsets %in% c("first", "all")) {
      if(offsets == "first") {
        foo = list(first_IFD)
        names(foo) = "first"
        attr(foo, "checksum") <- checksumXIF(fileName)
        attr(foo, "fileName_image") <- fileName
        attr(foo, "class") <- c("IFC_ifd_list", "IFC_first_ifd")
        return(foo)
      } else {
        offsets = suppressMessages(getOffsets(fileName, display_progress = display_progress, verbose = verbose, fast = TRUE))
      }
    }
  }
  
  K = c("IFC_ifd_list", "IFC_partial_ifd")
  if(!("IFC_offset" %in% class(offsets))) stop("'offsets' should be either 'all', 'first' or an object of class `IFC_offset`")
  L = length(offsets)
  if(L > 0) {
    if(L == obj_number*2) K = c("IFC_ifd_list", "IFC_full_ifd")
    VER = ifelse(verbose & (verbosity==2), TRUE, FALSE)
    if(display_progress) { 
      pb = newPB(session = dots$session, min = 0, max = L, initial = 0, style = 3)
      on.exit(endPB(pb), add = TRUE)
      ans = lapply(1:L, FUN=function(i_off) {
        setPB(pb, value = i_off, title = title_progress, label = "extracting IFDs")
        return(cpp_getTAGS(fname = fileName, offset = offsets[i_off], trunc_bytes = trunc_bytes, force_trunc = force_trunc, verbose = VER))
      }) 
    } else {
      ans = lapply(1:L, FUN=function(i_off) {
        return(cpp_getTAGS(fname = fileName, offset = offsets[i_off], trunc_bytes = trunc_bytes, force_trunc = force_trunc, verbose = VER))
      }) 
    }
    names(ans) = names(offsets)
  } else {
    ans = list()
    K = c("IFC_ifd_list", "IFC_empty_ifd")
  }
  attr(ans, "checksum") <- attr(offsets, "checksum")
  attr(ans, "fileName_image") <- fileName
  attr(ans, "class") <- K
  return(ans)
}
