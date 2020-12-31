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

#' @title IFC_pops XML Conversion
#' @description 
#' Helper to convert populations (`IFC_pops` object) to XML nodes.
#' @param pops an `IFC_pops` object.
#' @param verbose whether to display message about current action. Default is FALSE.
#' @param display_progress whether to display a progress bar. Default is TRUE.
#' @param title_progress character string, giving the title of the progress bar. Default is "".
#' @param ... other arguments to be passed.
#' @return a xml_node.
#' @keywords internal
toXML2_pops = function(pops, verbose = FALSE, display_progress = TRUE, title_progress = "", ...) {
  dots = list(...)
  assert(verbose, alw = c(TRUE, FALSE))
  if(verbose) message("creating pops node")
  assert(pops, cla = "IFC_pops")
  if(length(pops)==0) return(xml_new_node(name = "Pops", text = ""))
  assert(display_progress, alw = c(TRUE, FALSE))
  assert(title_progress, len = 1, typ = "character")
  tmp_style = c(20, 4, 3, 1, 5, 0, 2, 18, 15, 17)
  names(tmp_style)=c("Simple Dot","Cross","Plus","Empty Circle","Empty Diamond","Empty Square","Empty Triangle","Solid Diamond","Solid Square","Solid Triangle")
  L = length(pops)
  if(display_progress) {
    pb = newPB(session = dots$session, min = 0, max = L, initial = 0, style = 3)
    on.exit(endPB(pb))
    pops_nodes = xml_new_node(name = "Pops", .children = lapply(1:L, FUN=function(i_pop) {
      setPB(pb, value = i_pop, title = title_progress, label = "converting pops (xml)")
      pop = pops[[i_pop]]
      if(pop$color=="Cyan4") pop$color <- "Teal"
      if(pop$lightModeColor=="Cyan4") pop$lightModeColor <- "Teal"
      if(pop$color=="Green4") pop$color <- "Green"
      if(pop$lightModeColor=="Green4") pop$lightModeColor <- "Green"
      if(pop$color=="Chartreuse") pop$color <- "Lime"
      if(pop$lightModeColor=="Chartreuse") pop$lightModeColor <- "Lime"
      pop$style <- names(which(pop$style == tmp_style))[1]
      switch(pop$type,
             "B" = {
               xml_new_node(name = "Pop", attrs = pop[c("name", "type", "base", "color", "lightModeColor", "style")])
             },
             "C" = {
               xml_new_node(name = "Pop", attrs = pop[c("name", "type", "base", "color", "lightModeColor", "style", "definition")])
             },
             "G" = {
               if(length(pop$fy)==0) {
                 xml_new_node(name = "Pop", attrs = pop[c("name", "type", "base", "color", "lightModeColor", "style", "region", "fx")])
               } else { 
                 xml_new_node(name = "Pop", attrs = pop[c("name", "type", "base", "color", "lightModeColor", "style", "region", "fx", "fy")])
               }
             },
             "T" = {
               xml_new_node(name = "Pop", attrs = pop[c("name", "type", "base", "color", "lightModeColor", "style")],
                            .children = lapply(num_to_string(which(pop$obj)-1), FUN=function(o) {
                              paste0('ob O="', o, '"')
                            }))
             })
    }))
  } else {
    pops_nodes = xml_new_node(name = "Pops", .children = lapply(1:L, FUN=function(i_pop) {
      pop = pops[[i_pop]]
      if(pop$color=="Cyan4") pop$color <- "Teal"
      if(pop$lightModeColor=="Cyan4") pop$lightModeColor <- "Teal"
      if(pop$color=="Green4") pop$color <- "Green"
      if(pop$lightModeColor=="Green4") pop$lightModeColor <- "Green"
      if(pop$color=="Chartreuse") pop$color <- "Lime"
      if(pop$lightModeColor=="Chartreuse") pop$lightModeColor <- "Lime"
      pop$style <- names(which(pop$style == tmp_style))[1]
      switch(pop$type,
             "B" = {
               xml_new_node(name = "Pop", attrs = pop[c("name", "type", "base", "color", "lightModeColor", "style")])
             },
             "C" = {
               xml_new_node(name = "Pop", attrs = pop[c("name", "type", "base", "color", "lightModeColor", "style", "definition")])
             },
             "G" = {
               if(length(pop$fy)==0) {
                 xml_new_node(name = "Pop", attrs = pop[c("name", "type", "base", "color", "lightModeColor", "style", "region", "fx")])
               } else {
                 xml_new_node(name = "Pop", attrs = pop[c("name", "type", "base", "color", "lightModeColor", "style", "region", "fx", "fy")])
               }
             },
             "T" = {
               xml_new_node(name = "Pop", attrs = pop[c("name", "type", "base", "color", "lightModeColor", "style")],
                            .children = lapply(num_to_string(which(pop$obj)-1), FUN=function(o) {
                              paste0('ob O="', o, '"')
                            }))
             })
    }))
  }
  return(pops_nodes)
}
