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

#' @title IFC_graphs XML Conversion
#' @description 
#' Helper to convert graphs (`IFC_graphs` object) to XML nodes.
#' @param graphs an `IFC_graphs` object.
#' @param verbose whether to display message about current action. Default is FALSE.
#' @return a xml_node.
#' @keywords internal
toXML2_graphs = function(graphs, verbose = FALSE) {
  assert(verbose, alw = c(TRUE, FALSE))
  if(verbose) message("creating graphs node")
  assert(graphs, cla = "IFC_graphs")
  if(length(graphs)==0) return(xml_new_node(name = "Displays", text = ""))
  graphs = lapply(graphs, FUN=function(i) {
    i$GraphRegion = lapply(i$GraphRegion, FUN = function(g) g[!grepl("def", names(g))]) # it is mandatory to remove def
    if(typeof(i) %in% c("integer","double")) {
      return(num_to_string(i))
    } else {
      return(i)
    }
  })
  xml_new_node(name = "Displays", attrs = list(count = num_to_string(length(graphs)), layout="1", entriesPerRow="12", lightDarkMode="1"),
             .children = lapply(graphs, FUN=function(i_graph) {
               xml_new_node(name = "Graph", attrs = i_graph[!grepl("Legend|BasePop|GraphRegion|ShownPop", names(i_graph)) & sapply(i_graph, FUN=function(ele) length(ele)!=0)],
                          .children = c(lapply(i_graph[["Legend"]], FUN=function(i) xml_new_node(name = "Legend", attrs = i)),
                                        lapply(i_graph[["BasePop"]], FUN=function(i) xml_new_node(name = "BasePop", attrs = i)),
                                        lapply(i_graph[["GraphRegion"]], FUN=function(i) xml_new_node(name = "GraphRegion", attrs = i)),
                                        lapply(i_graph[["ShownPop"]], FUN=function(i) xml_new_node(name = "ShownPop", attrs = i))))
  }))
}
