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

#' @title Remove Features from an IFC_data Object
#' @description
#' Removes regions from an already existing `IFC_data` object.
#' @param obj an `IFC_data` object extracted by ExtractFromDAF(extract_features = TRUE) or ExtractFromXIF(extract_features = TRUE).
#' @param features a character vector of features names to remove within 'obj'. Note that "Object Number" is not allowed and will be excluded from 'features' if present.
#' @param list_only whether to return a list of elements that will be impacted by the removal. Default is TRUE.
#' If FALSE then modified object will be returned
#' @param ... Other arguments to be passed.
#' @return an `IFC_data` object or a list of elements impacted by removal depending on 'list_only' parameter.
#' @export
data_rm_features <- function(obj, features, list_only = TRUE, ...) {
  dots = list(...)
  assert(obj, cla = "IFC_data")
  assert(list_only, len = 1, alw = c(TRUE,FALSE))
  assert(features, typ = "character")
  to_remove_features = features
  if(length(obj$features) == 0) {
    warning("'obj' contains no feature", immediate. = TRUE, call. = FALSE)
    if(list_only) {
      return(list(masks = character(),
                  features = character(),
                  regions = character(),
                  pops = character(),
                  graphs = integer()))
    } else {
      return(obj)
    }
  }
  
  # forbids removal of "Objet Number"
  if(any("Object Number" %in% to_remove_features)) warning("\"Object Number\" is mandatory and can not be removed", immediate. = TRUE, call. = FALSE)
  to_remove_features = setdiff(to_remove_features, c("Object Number"))
  
  # removes duplicated inputs
  tmp = duplicated(to_remove_features)
  if(any(tmp)) {
    warning(paste0("duplicated 'features' automatically removed:\n", paste0(paste0("\t- ", to_remove_features[tmp]), collapse = "\n")), immediate. = TRUE, call. = FALSE)
    to_remove_features = to_remove_features[!tmp]
  }
  
  # removes features not in obj
  tmp = to_remove_features %in% names(obj$features)
  if(any(!tmp)) {
    warning(paste0("some 'features' are not in 'obj$features' and can't be removed:\n", paste0(paste0("\t- ", to_remove_features[!tmp]), collapse = "\n")), immediate. = TRUE, call. = FALSE)
    to_remove_features = to_remove_features[tmp]
  }
  if(length(to_remove_features) == 0) {
    warning("no feature to remove in 'obj'", immediate. = TRUE, call. = FALSE)
    if(list_only) {
      return(list(masks = character(),
                  features = character(),
                  regions = character(),
                  pops = character(),
                  graphs = integer()))
    } else {
      return(obj)
    }
  }
  
  # search for features depending on input features
  L = length(to_remove_features)
  LL = L - 1
  while(L != LL) {
    L = length(to_remove_features)
    for(i in 1:length(obj$features_def)) {
      if((obj$features_def[[i]]$type == "combined") && any(to_remove_pops %in% strsplit(obj$features_def[[i]]$def, split = "|", fixed = TRUE)[[1]])) {
        to_remove_features = c(to_remove_features, obj$features_def[[i]]$name)
      }
    }
    to_remove_features = unique(to_remove_features)
    LL = length(to_remove_features)
  }
  
  # search regions that depend on input features
  to_remove_regions = character()
  for(i in 1:length(obj$pops)) {
    if(obj$pops[[i]]$type == "G" && any(to_remove_features %in% unlist(obj$pops[[i]][c("fx", "fy")]))) to_remove_regions = c(to_remove_regions, obj$pops[[i]]$region)
  }
  to_remove_regions = unique(to_remove_regions)
  
  # search pops that depend on input features
  to_remove_pops = character()
  for(i in 1:length(obj$pops)) {
    if((obj$pops[[i]]$type == "G") && any(to_remove_regions %in% obj$pops[[i]]$region)) to_remove_pops = c(to_remove_pops, obj$pops[[i]]$name)
  }
  to_remove_pops = unique(to_remove_pops)
  
  # search for pops that depend on previous pops
  for(i in 1:length(obj$pops)) {
    if(any(to_remove_pops %in% c(obj$pops[[i]]$base, obj$pops[[i]]$names))) to_remove_pops = c(to_remove_pops, obj$pops[[i]]$name)
  }
  to_remove_pops = unique(to_remove_pops)
  
  # search graphs that depend on input pops
  to_remove_graphs = integer()
  if(length(obj$graphs) > 0) for(i in 1:length(obj$graphs)) {
    g = obj$graphs[[i]]
    base = sapply(g$BasePop, FUN = function(p) p$name)
    shown = sapply(g$ShownPop, FUN = function(p) p$name)
    region = sapply(g$GraphRegion, FUN = function(r) r$name)
    region_def = sapply(g$GraphRegion, FUN = function(r) r$def)
    if(any(c(to_remove_pops %in% c(base, shown, region_def), to_remove_features %in% unlist(g[c("f1","f2")])))) {
      to_remove_graphs = c(to_remove_graphs, i)
    }
  }
  to_remove_graphs = unique(to_remove_graphs)
  
  # create list
  if(list_only) {
    return(list(masks = character(),
                features = to_remove_features,
                regions = to_remove_regions,
                pops = to_remove_pops,
                graphs = to_remove_graphs))
  }
  
  # remove regions and their dep
  if(length(to_remove_features) != 0) {
    tmp = names(obj$features) %in% to_remove_features
    if(any(!tmp)) {
      if(sum(!tmp) != 1) {
        obj$features = structure(obj$features[, !tmp], class = class(obj$features))
      } else {
        obj$features = structure(as.data.frame(obj$features[, !tmp]), class = class(obj$features), stringsAsFactors = FALSE)
      }
      names(obj$features) = names(obj$features_def[!tmp])
      obj$features_def = structure(obj$features_def[!tmp], class = class(obj$features_def))
    } else {
      obj$features = structure(data.frame(), class = class(obj$features))
      obj$features_def = structure(list(), class = class(obj$features_def))
    }
  }
  
  if(length(to_remove_regions) != 0) {
    tmp = names(obj$regions) %in% to_remove_regions
    if(any(!tmp)) {
      obj$regions = structure(obj$regions[!tmp], class = class(obj$regions))
    } else {
      obj$regions = structure(list(), class = class(obj$regions))
    }
  }
  pops_back = obj$pops
  obj$pops = list()
  obj = data_add_pops(obj, pops = pops_back[!(names(pops_back) %in% to_remove_pops)], ...)
  if(length(to_remove_graphs) != 0) return(adjustGraph(obj = obj, selection = to_remove_graphs))
  return(obj)
}
