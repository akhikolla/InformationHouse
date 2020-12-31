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

#' @title IFC_images XML Conversion
#' @description 
#' Helper to convert images (`IFC_images` object) to XML nodes.
#' @param images an `IFC_images` object.
#' @param verbose whether to display message about current action. Default is FALSE.
#' @param title_progress character string, giving the title of the progress bar. Default is "".
#' @return a xml_node.
#' @keywords internal
toXML2_images = function(images, verbose = FALSE, display_progress = TRUE, title_progress = "") {
  assert(verbose, alw = c(TRUE, FALSE))
  if(verbose) message("creating images node")
  assert(images, cla = "IFC_images")
  bgm = grep("^bgmean", names(images))
  bgs = grep("^bgstd", names(images))
  satc = grep("^satcount", names(images))
  satp = grep("^satpercent", names(images))
  lapply(1:nrow(images), FUN=function(i) {
    xml_new_node(name = "SO",
               attrs = c("id" = num_to_string(images[i, "id"]),
                         "imgIFD" = num_to_string(images[i, "imgIFD"]),
                         "mskIFD" = num_to_string(images[i, "mskIFD"]),
                         "spIFD" = num_to_string(images[i, "spIFD"]),
                         "w" = num_to_string(images[i, "w"]),
                         "l" = num_to_string(images[i, "l"]),
                         "fs" = num_to_string(images[i, "fs"]),
                         "cl" = num_to_string(images[i, "cl"]),
                         "ct" = num_to_string(images[i, "ct"]),
                         "objCenterX" = num_to_string(images[i, "objCenterX"]),
                         "objCenterY" = num_to_string(images[i, "objCenterY"]),
                         "bgmean" = paste0(num_to_string(unlist(images[i, bgm])), collapse = "|"),
                         "bgstd" = paste0(num_to_string(unlist(images[i, bgs])), collapse = "|"),
                         "satcount" = paste0(num_to_string(unlist(images[i, satc])), collapse = "|"),
                         "satpercent" = paste0(num_to_string(unlist(images[i, satp])), collapse = "|")))
    
  })
}
