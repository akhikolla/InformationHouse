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

#' @title IFC Files Generic Writer
#' @description
#' Writes IFC data to DAF and subsets or merges RIF/CIF Files.
#' @param fileName path to file.
#' @param ... arguments to pass to \code{\link{ExportToDAF}} or \code{\link{ExportToXIF}}.
#' @details If 'fileName' is a DAF file \code{\link{ExportToDAF}} will be used to write file whereas if it is a RIF or CIF file \code{\link{writeIFC}} will use \code{\link{ExportToXIF}}.
#' @return it invisible returns the path of exported file.
#' @examples
#' if(requireNamespace("IFCdata", quietly = TRUE)) {
#'   tmp <- tempdir(check = TRUE)
#'   ## use a daf file
#'   file_daf <- system.file("extdata", "example.daf", package = "IFCdata")
#'   ## create a tagged population named test with 1st object
#'   pop <- buildPopulation(name = "test", type = "T", obj = 0)
#'   writeIFC(file_daf, write_to = paste0(tmp, "\\test_write.daf"),
#'            overwrite = TRUE, pops = list(pop))
#'   ## use a rif file, but you can also use a cif
#'   file_rif <- system.file("extdata", "example.rif", package = "IFCdata")
#'   writeIFC(fileName = file_rif, write_to = paste0(tmp, "\\test_write.rif"), 
#'            overwrite = TRUE, objects = 0)
#' } else {
#'   message(sprintf('Please run `install.packages("IFCdata", repos = "%s", type = "source")` %s',
#'                   'https://gitdemont.github.io/IFCdata/',
#'                   'to install extra files required to run this example.'))
#' }
#' @export
writeIFC <- function(fileName, ...) {
  dots=list(...)
  file_extension = getFileExt(fileName)
  assert(file_extension, alw = c("daf", "rif", "cif"))
  if((length(file_extension) == 1) && (file_extension == "daf")) return(ExportToDAF(fileName = fileName, ...))
  return(ExportToXIF(fileName = fileName, ...))
}
