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

#' @title IFC_regions XML Conversion
#' @description 
#' Helper to convert regions (`IFC_regions` object) to XML nodes.
#' @param regions an `IFC_regions` object.
#' @param verbose whether to display message about current action. Default is FALSE.
#' @return a xml_node.
#' @keywords internal
toXML2_regions = function(regions, verbose = FALSE) {
  assert(verbose, alw = c(TRUE, FALSE))
  if(verbose) message("creating regions node")
  assert(regions, cla = "IFC_regions")
  if(length(regions)==0) return(xml_new_node(name = "Regions", text = ""))
  xml_new_node(name = "Regions", .children = lapply(regions, FUN=function(i_reg) {
    if(i_reg$color=="Cyan4") i_reg$color <- "Teal"
    if(i_reg$lightcolor=="Cyan4") i_reg$lightcolor <- "Teal"
    if(i_reg$color=="Green4") i_reg$color <- "Green"
    if(i_reg$lightcolor=="Green4") i_reg$lightcolor <- "Green"
    if(i_reg$color=="Chartreuse") i_reg$color <- "Lime"
    if(i_reg$lightcolor=="Chartreuse") i_reg$lightcolor <- "Lime"
    xml_new_node(name = "Region",
               attrs = i_reg[!grepl("^x$|^y$", names(i_reg))],
               .children = lapply(1:length(i_reg[["x"]]), FUN = function(i_coord) {
                 xml_new_node(name = "axy", attrs = list(x = num_to_string(i_reg[["x"]][i_coord]), y = num_to_string(i_reg[["y"]][i_coord])))
               }))
  }))
}
