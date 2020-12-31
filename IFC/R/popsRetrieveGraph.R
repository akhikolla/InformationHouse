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

#' @title Graph Retrieval from Graphical IFC_pops
#' @description
#' Retrieves the graph a graphical population originate from
#' @param obj an `IFC_data` object extracted by ExtractFromDAF(extract_features = TRUE) or ExtractFromXIF(extract_features = TRUE).
#' @param pops names of graphical populations present in 'obj'. Note that they should be siblings.
#' @param vis2D when original graph is not an histogram, whether to display it as "scatter" or "density". Default is "density".
#' @param all_siblings whether to add all 'pop' siblings in the graph. Default is FALSE.
#' @return a list of parameters needed to build an IFC graph.
#' @keywords internal
popsRetrieveGraph = function(obj, pops, vis2D = "density", all_siblings = FALSE) {
  if(missing(obj)) stop("'obj' can't be missing")
  if(!("IFC_data"%in%class(obj))) stop("'obj' is not of class `IFC_data`")
  if(missing(pops)) stop("'pops' can't be missing")
  pops = unique(pops); pops = as.character(pops); assert(pops, alw = names(obj$pops))
  if(!all(sapply(obj$pops[pops], FUN=function(p) p$type=="G"))) stop("'pops' should be of type graphical")
  # vis2D = as.character(vis2D); assert(vis2D, len=1, alw=c("scatter","density"))
  all_siblings = as.logical(all_siblings); assert(all_siblings, len = 1, alw =c(TRUE, FALSE))
  
  siblings = popsGetSiblings(obj, pops)
  if(length(siblings) == 0) stop("'pops' should be siblings")
  
  # initializes variables
  if(all_siblings) {
    pops = siblings
  } else {
    pops = pops
  }
  sib1 = popsGetSiblings1(obj, pops)
  sib2 = popsGetSiblings2(obj, pops)
  parent1 = lapply(obj$pops[names(sib1)], FUN = function(p) p$base)
  parent1 = unique(unlist(parent1))
  parent2 = lapply(obj$pops[names(sib2)], FUN = function(p) p$base)
  parent2 = unique(unlist(parent2))

  P = obj$pops[pops]
  SUB = apply(do.call("rbind", lapply(obj$pops[unique(c(parent1, parent2))], FUN = function(p) p$obj)), 2, any)
  R = sapply(P, simplify = F, FUN=function(p) obj$regions[[p$region]])
  foo = list()
  # start rebuilding original graph
  foo$f1 = P[[1]]$fx
  foo$xlogrange = R[[1]]$xlogrange
  foo$ShownPop = list()
  foo$title = paste0(unique(c(parent1, parent2)), collapse = ", ")
  if(length(P[[1]]$fy) == 0) {
    xran = range(c(obj$features[SUB, foo$f1], unlist(lapply(R, FUN=function(r) c(r$x, r$cx)))), na.rm = TRUE)
    if(foo$xlogrange == "P") {
      xran = xran + diff(xran) * c(-0.07,0.07)
    } else {
      xran = smoothLinLog(xran, hyper = as.numeric(foo$xlogrange))
      xran = xran + diff(xran) * c(-0.07,0.07)
    }
    foo$xmin = xran[1]
    foo$xmax = xran[2]
    foo$type = "histogram"
    foo$bincount = 0
    foo$freq = "T"
    br = do.breaks(xran, 520)
    yran = c(0,max(sapply(obj$pops[unique(c(parent1, parent2))], FUN=function(p) {
      x = obj$features[p$obj, foo$f1]
      if(foo$xlogrange != "P") x = smoothLinLog(x, hyper = as.numeric(foo$xlogrange))
      get_ylim(x=x, type="percent", br=br) * 1.07
    })))
    if(yran[1] == yran[2]) yran = yran[1] + c(0,0.07)
  } else {
    xran = range(c(obj$features[SUB, foo$f1], unlist(lapply(R, FUN=function(r) c(r$x, r$cx)))), na.rm = TRUE)
    if(foo$xlogrange == "P") {
      xran = xran + diff(xran) * c(-0.07,0.07)
    } else {
      xran = smoothLinLog(xran, hyper = as.numeric(foo$xlogrange))
      xran = xran + diff(xran) * c(-0.07,0.07)
    }
    foo$f2 = P[[1]]$fy
    foo$ylogrange = R[[1]]$ylogrange
    yran = range(c(obj$features[SUB, foo$f2], unlist(lapply(R, FUN=function(r) c(r$y,r$cy)))), na.rm = TRUE)
    if(foo$ylogrange == "P") {
      yran = yran + diff(yran) * c(-0.07,0.07)
    } else {
      yran = smoothLinLog(yran, hyper = as.numeric(foo$ylogrange))
      yran = yran + diff(yran) * c(-0.07,0.07)
      yran = inv_smoothLinLog(yran, hyper = as.numeric(foo$ylogrange))
    }
    foo$type = vis2D
  }
  if(foo$xlogrange != "P") xran = inv_smoothLinLog(xran, hyper = as.numeric(foo$xlogrange))
  foo$xmin = xran[1]
  foo$xmax = xran[2]
  foo$ymin = yran[1]
  foo$ymax = yran[2]
  foo$BasePop = lapply(unique(c(parent1, parent2)), FUN = function(p) list(name = p, linestyle = "Solid", fill = "true"))
  foo$GraphRegion = list()
  if(length(R) > 0) foo$GraphRegion = list(list("name" = R[[1]]$label, def = c(R[[1]]$def, names(R)[1])))
  if(length(R) > 1) for(i_reg in 2:length(R)) {
    defined = sapply(foo$GraphRegion, FUN = function(r) r$name) %in% R[[i_reg]]$label
    if(any(defined)) {
      foo$GraphRegion[[defined]] = list("name" = R[[i_reg]]$label, def = c(foo$GraphRegion[[defined]]$def, names(R)[i_reg]))
    } else {
      foo$GraphRegion = c(foo$GraphRegion, list(list("name" = R[[i_reg]]$label, def = names(R)[i_reg])))
    }
  }
  return(foo)
}
