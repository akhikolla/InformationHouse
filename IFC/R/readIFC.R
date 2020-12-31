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

#' @title IFC Files Generic Reader
#' @description
#' Reads IFC data from IFC files no matter if they are DAF, RIF or CIF.
#' @param fileName path to file.
#' @param ... arguments to pass to \code{\link{ExtractFromDAF}} or \code{\link{ExtractFromXIF}}.
#' @details If input 'fileName' is a DAF file \code{\link{ExtractFromDAF}} will be used to read the file whereas if it is a CIF or RIF file \code{\link{readIFC}} will use \code{\link{ExtractFromXIF}}.
#' @return an object of class `IFC_data`.
#' @examples
#' if(requireNamespace("IFCdata", quietly = TRUE)) {
#'   ## use a rif file, but you can also read daf or cif
#'   file_rif <- system.file("extdata", "example.rif", package = "IFCdata")
#'   rif <- readIFC(fileName = file_rif)
#' } else {
#'   message(sprintf('Please run `install.packages("IFCdata", repos = "%s", type = "source")` %s',
#'                   'https://gitdemont.github.io/IFCdata/',
#'                   'to install extra files required to run this example.'))
#' }
#' @export
readIFC <- function(fileName, ...) {
  dots=list(...)
  file_extension = getFileExt(fileName)
  assert(file_extension, len = 1, alw = c("daf", "rif", "cif"))
  if(file_extension == "daf") return(ExtractFromDAF(fileName = fileName, ...))
  return(ExtractFromXIF(fileName = fileName, ...))
}
