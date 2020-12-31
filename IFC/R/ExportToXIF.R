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

#' @title RIF/CIF File Writer
#' @description
#' Subsets or merges RIF or CIF files.
#' @param fileName path(s) of file(s) to subset or merge.
#' If multiple files are provided they will be merged. Otherwise, if only one file is input it will be subsetted.
#' All files have to be either '.rif' or '.cif' files.
#' All files should have same channels.
#' @param write_to pattern used to export file.
#' Placeholders, like "\%d/\%s_fromR.\%e", will be substituted:\cr
#' -\%d: with full path directory of first element of 'fileName'\cr
#' -\%p: with first parent directory of first element of 'fileName'\cr
#' -\%e: with extension of first element of 'fileName' (without leading .)\cr
#' -\%s: with shortname from of first element of 'fileName' (i.e. basename without extension).\cr
#' Exported file extension will be deduced from this pattern. It has to be the same as 'fileName', i.e. .cif or .rif.
#' @param objects integer vector, IDEAS objects ids numbers to use. If missing, the default, all objects will be used. Only apply for subsetting.
#' @param offsets object of class `IFC_offset`. If missing, the default, offsets will be extracted from 'fileName'.\cr
#' This param is not mandatory but it may allow to save time for repeated XIF export on same file. Only apply for subsetting.
#' @param fast whether to fast extract 'objects' or not. Default is TRUE.
#' Meaning that 'objects' will be extracting expecting that objects are stored in ascending order.\cr
#' Note that a warning will be sent if an 'object' is found at an unexpected order.
#' In such a case you may need to rerun function with 'fast' = FALSE.
#' If set to FALSE, all object_ids will be scanned from 'fileName' to ensure extraction of desired 'objects'.\cr
#' IMPORTANT: whatever this argument is, features are extracted assuming an ascending order of storage in file.\cr
#' Only apply for subsetting.
#' @param extract_features whether to try to extract features. Default is FALSE.
#' IMPORTANT: it is not clear if how features are stored and which objects they rely to when input file is already a merge or a subset.
#' For this reason it should be carefully checked.
#' Note that features extraction is not implemented for merging.
#' @param endianness the endian-ness ("big" or "little") of the target system for the file. Default is .Platform$endian.\cr
#' Endianness describes the bytes order of data stored within the files. This parameter may not be modified.
#' @param verbose whether to display information (use for debugging purpose). Default is FALSE.
#' @param verbosity quantity of information displayed when verbose is TRUE; 1: normal, 2: rich. Default is 1.
#' @param overwrite whether to overwrite file or not. Default is FALSE.\cr
#' Note that if TRUE, it will overwrite exported file if path(s) of file(s) in 'fileName' and deduced from 'write_to' arguments are different.
#' Otherwise, you will get an error saying that overwritting source file is not allowed.\cr
#' Note also that an original file, i.e. generated by IDEAS(R) or INSPIRE(R), will never be overwritten.
#' Otherwise, you will get an error saying that overwritting original file is not allowed.\cr
#' @param display_progress whether to display a progress bar. Default is TRUE.
#' @param ... other arguments to be passed.
#' @details when 'extract_features' is set TRUE, only features stored in binary format will be extracted if found.\cr
#' If the input 'fileName' is a merged of several files then features will be extracted from these files.\cr
#' If these files can't be found, Warning(s) will be thrown and input 'fileName' will be extracted without features values.
#' @examples
#' if(requireNamespace("IFCdata", quietly = TRUE)) {
#'   tmp <- tempdir(check = TRUE)
#'   ## use a cif file, but you can also subset rif
#'   file_cif <- system.file("extdata", "example.cif", package = "IFCdata")
#'   ## subset objects 0,1 and 4 from file
#'   exported <- ExportToXIF(fileName = file_cif, write_to = paste0(tmp, "\\test.cif"),
#'                           overwrite = TRUE, objects = c(0,1,4))
#' } else {
#'   message(sprintf('Please run `install.packages("IFCdata", repos = "%s", type = "source")` %s',
#'                   'https://gitdemont.github.io/IFCdata/',
#'                   'to install extra files required to run this example.'))
#' }
#' @return It invisibly returns full path of exported file.
#' @export
ExportToXIF <- function(fileName, write_to, 
                        objects, offsets, fast = TRUE, 
                        extract_features = FALSE, endianness = .Platform$endian, verbose = FALSE, verbosity = 1, 
                        overwrite = FALSE, display_progress = TRUE, ...) {
  dots = list(...)
  # check madatory param
  if(missing(fileName)) stop("'fileName' can't be missing")
  if(missing(write_to)) stop("'write_to' can't be missing")
  tmp = duplicated(fileName)
  if(any(tmp)) {
    warning(paste0("duplicated files have been removed from 'fileName': ","\n-", paste0(fileName[tmp],collapse="\n-")))
    fileName = fileName[!tmp]
  }
  args = list(fileName = fileName, write_to = write_to,
              extract_features = extract_features, endianness = endianness, verbose = verbose, verbosity = verbosity, 
              overwrite = overwrite, display_progress = display_progress)
  args = c(args, dots)
  if(length(fileName) > 1) {
    return(do.call(what = mergeXIF, args = args))
  } else {
    if(!missing(offsets)) args = c(args, list(offsets = offsets))
    if(!missing(objects)) args = c(args, list(objects = objects))
    return(do.call(what = subsetXIF, args = args))
  }
}
