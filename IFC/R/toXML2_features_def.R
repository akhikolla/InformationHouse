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

#' @title IFC_features_def XML Conversion
#' @description 
#' Helper to convert features definition (`IFC_features_def` object) to XML nodes.
#' @param features an `IFC_features_def` object.
#' @param verbose whether to display message about current action. Default is FALSE.
#' @return a XML::xmlNode.
#' @keywords internal
toXML2_features_def =function(features_def, verbose = verbose) {
  assert(verbose, alw = c(TRUE, FALSE))
  if(verbose) message("creating features definition node")
  assert(features_def, cla = "IFC_features_def")
  if(length(features_def)==0) return(xml_new_node(name = "DefinedFeatures", text = ""))
  xml_new_node(name = "DefinedFeatures", .children = lapply(features_def, FUN=function(i_def) {
    xml_new_node(name = "UDF", attrs = i_def)
  }))
}
