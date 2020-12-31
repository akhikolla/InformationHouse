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

#' @title IFC Files Checksum
#' @description
#' This function returns RIF/CIF checksum.
#' Checksum is the sum of img IFDs (Image Field Directory) offsets of objects 0, 1, 2, 3 and 4.
#' @param fileName path to file.
#' @param ... arguments to pass to \code{\link{checksumDAF}} or \code{\link{checksumXIF}}.
#' @details if fileName is a DAF file, then CIF checksum is computed from images values found in DAF.
#' @return an integer corresponding to IFC file checksum.
#' @export
checksumIFC <- function(fileName, ...) {
  dots=list(...)
  file_extension = getFileExt(fileName)
  assert(file_extension, len = 1, alw = c("daf", "rif", "cif"))
  if(file_extension == "daf") return(checksumDAF(fileName = fileName, ...))
  return(checksumXIF(fileName = fileName, ...))
}
