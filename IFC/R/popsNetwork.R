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

#' @title IFC_pops Network Display
#' @description
#' Builds and displays populations network.
#' @param obj an `IFC_data` object extracted with features extracted.
#' @param hierarchical whether to display network using a hierarchical layout or not. Default is TRUE.
#' @param color_mode Whether to extract colors from 'obj' in "white" or "black" mode. Default is "white".
#' @param highlight population to permanently highlight. If found in 'obj', this population will be displayed with its color. Default is NULL.
#' @param seed If you provide a seed manually, the layout will be the same every time. Default is NULL.
#' @param direction The direction of the hierarchical layout. Default is 'LR'.\cr
#' The available options are: 'UD', 'DU', 'LR', 'RL'. To simplify: up-down, down-up, left-right, right-left.
#' @param weighted whether to scale population's node size according to count. Default is TRUE.
#' @param ... other argument to be passed.
#' @examples
#' if(requireNamespace("IFCdata", quietly = TRUE)) {
#'   ## use a daf file
#'   file_daf <- system.file("extdata", "example.daf", package = "IFCdata")
#'   daf <- ExtractFromDAF(fileName = file_daf)
#'   popsNetwork(obj = daf)
#' } else {
#'   message(sprintf('Please run `install.packages("IFCdata", repos = "%s", type = "source")` %s',
#'                   'https://gitdemont.github.io/IFCdata/',
#'                   'to install extra files required to run this example.'))
#' }
#' @return a \pkg{visNetwork} object.
#' @export
popsNetwork = function(obj, hierarchical=TRUE, color_mode="white", highlight=NULL, seed=NULL, direction="LR", weighted=TRUE, ...) {
  dots = list(...)
  if(missing(obj)) stop("'obj' can't be missing")
  if(!("IFC_data"%in%class(obj))) stop("'obj' is not of class `IFC_data`")
  if(length(obj$pops)==0) stop("please use argument 'extract_features' = TRUE with ExtractFromDAF() or ExtractFromXIF() and ensure that features were correctly extracted")
  if(!is.null(seed)) {seed=na.omit(as.integer(seed)); assert(seed, len=1, typ="integer")}
  hierarchical = as.logical(hierarchical); assert(hierarchical, len=1, alw=c(TRUE,FALSE))
  weighted = as.logical(weighted); assert(weighted, len=1, alw=c(TRUE,FALSE))
  assert(color_mode, len=1, alw=c("black","white"))
  assert(direction, len=1, alw=c("UD", "DU", "LR", "RL"))
  color_mode=which(c("black","white")%in%color_mode)
  # determine levels for hierarchical layout
  lev = popsGetLevels(obj$pops)
  size = table(lev)
  size = c(length(size), max(size))
  width = NULL
  height = NULL
  if(hierarchical && length(dots$session) != 0) {
    width = size[1] * 75
    height = size[2] * 75
  }
  # create nodes
  nodes = as.data.frame(t(sapply(obj$pops, FUN=function(p) {
    typ = c("Base","Boolean","Graphical","Tagged")[c("B","C","G","T")%in%p$type]
    col = rgb(t(checkColor(p[c("color","lightModeColor")][[color_mode]])), maxColorValue = 255)
    val = ifelse(weighted, sum(p$obj), 1)
    tit = paste0("<em>", p$name, "</em><br>", typ, ", n: <b>", val, "</b>, %All: <b>", round(100*sum(p$obj)/length(p$obj),2), "</b>")
    if(typ=="Boolean") tit = paste0(tit, "<br>def: ", paste0(p$split,collapse=" "))
    if(typ=="Graphical") {
      tit = paste0(tit, ", %", p$base, ": <b>", round(100*sum(p$obj)/sum(obj$pops[[p$base]]$obj),2), "</b>")
      tit = paste0(tit, "<br><b>", obj$regions[[p$region]]$type, "</b><br><b>x</b>: ", p$fx)
      if(!is.null(p$fy)) tit = paste(tit, "<br><b>y</b>: ", p$fy)
      tit = paste0(tit, "<table> <tr> <th>x</th> <th>y</th> </tr>")
      tit = paste0(tit, paste("<tr>", apply(do.call(cbind, args = obj$regions[[p$region]][c("x","y")]), 1, FUN=function(r) paste("<td>",r,"</td>",collapse=" ")), "</tr>", collapse = "\n"), "</table>")
    }
    c("label"=p$name, "id"=p$name, "type"=typ, "title"=tit, "value"=val, 
      "color.background"=ifelse(p$name%in%highlight, col, "lightgrey"),
      "color.border"=ifelse(!p$name%in%highlight, col, "lightgrey"),
      "color.highlight.background"=col, 
      "color.highlight.border"=ifelse(!p$name%in%highlight, col, "lightgrey"),
      "color.hover.background"=col,
      "color.hover.border"=ifelse(!p$name%in%highlight, col, "lightgrey"))
  })), stringsAsFactors = FALSE)
  nodes$level = lev
  # create edges
  edges = lapply(obj$pops, FUN=function(p) {
    col = rgb(t(checkColor(p[c("color","lightModeColor")][[color_mode]])), maxColorValue = 255)
    if(p$type=="C") {
      from = p$names
    } else {
      from = p$base
    }
    return(data.frame(to=rep(p$name, length(from)), from=from, arrows="to", stringsAsFactors = FALSE))
  })
  visNetwork(nodes, do.call(rbind, edges[-1]), main = list(text="Populations network"), width = width, height = height) %>%
    visPhysics(enabled = TRUE, stabilization = FALSE) %>%
    visOptions(highlightNearest = list(enabled=TRUE, hover=TRUE, labelOnly=FALSE, degree=list(from=1,to=1), algorithm="hierarchical"), selectedBy="type") %>%
    visLayout(randomSeed = seed) %>%
    visInteraction(navigationButtons = TRUE) %>%
    visHierarchicalLayout(direction = direction, enabled = hierarchical) 
}
