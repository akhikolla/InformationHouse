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

#' @title Add Region to IFC_data Object
#' @description
#' Adds regions to an already existing `IFC_data` object.
#' @param obj an `IFC_data` object extracted by ExtractFromDAF(extract_features = TRUE) or ExtractFromXIF(extract_features = TRUE).
#' @param regions a list of region(s) to add to obj. Each element of this list will be coerced by \code{\link{buildRegion}}.
#' @details A warning will be thrown if a provided region is already existing in 'obj'.\cr
#' In such a case this region will not be added to 'obj'.\cr
#' If any input population is not well defined and can't be created then an error will occur.
#' @param ... Other arguments to be passed.
#' @examples
#' if(requireNamespace("IFCdata", quietly = TRUE)) {
#'   ## use a daf file
#'   file_daf <- system.file("extdata", "example.daf", package = "IFCdata")
#'   daf <- ExtractFromDAF(fileName = file_daf)
#'   ## copy 1st region found in daf
#'   reg <- daf$regions[[1]]
#'   if(length(reg) != 0) {
#'     reg_copy <- reg
#'     ## modify region label and x boundaries
#'     reg_copy$label <- paste0(reg_copy$label,"_copy")
#'     reg_copy$x <- reg_copy$x*0.9
#'     ## create new object with this new region
#'     dafnew <- data_add_regions(obj = daf, regions = list(reg_copy))
#'   }
#' } else {
#'   message(sprintf('Please run `install.packages("IFCdata", repos = "%s", type = "source")` %s',
#'                   'https://gitdemont.github.io/IFCdata/',
#'                   'to install extra files required to run this example.'))
#' }
#' @return an IFC_data object with regions added.
#' @export
data_add_regions <- function(obj, regions, ...) {
  assert(obj, cla = "IFC_data")
  
  # try to coerce regions inputs
  regions = lapply(regions, FUN=function(x) do.call(what=buildRegion, args=x))
  names(regions) = sapply(regions, FUN=function(x) x$label)
  
  # removes duplicated inputs
  tmp = duplicated(names(regions))
  if(any(tmp)) {
    warning(paste0("duplicated regions automatically removed: ", names(regions)[tmp]), immediate. = TRUE, call. = FALSE)
    regions = regions[!tmp]
  }
  
  # change colors to R compatible
  for(i in 1:length(regions)) {
    if(regions[[i]]$color=="Teal") regions[[i]]$color <- "Cyan4"
    if(regions[[i]]$color=="Green") regions[[i]]$color <- "Green4"
    if(regions[[i]]$color=="Lime") regions[[i]]$color <- "Chartreuse"
    if(regions[[i]]$lightcolor=="Teal") regions[[i]]$lightcolor <- "Cyan4"
    if(regions[[i]]$lightcolor=="Green") regions[[i]]$lightcolor <- "Green4"
    if(regions[[i]]$lightcolor=="Lime") regions[[i]]$lightcolor <- "Chartreuse"
  }
  
  # removes already defined regions
  exported_regions = sapply(regions, FUN=function(reg) {
    if(reg$label%in%names(obj$regions)) {
      warning(paste0(reg$label, "\nnot exported: trying to export an already defined region"), immediate. = TRUE, call. = FALSE)
      return(FALSE)
    }
    return(TRUE)
  })
  exported_regions = regions[exported_regions]
  names(exported_regions) = sapply(exported_regions, FUN=function(x) x$label)
  
  K = class(obj$regions)
  obj$regions = c(obj$regions, exported_regions)
  class(obj$regions) = K
  
  return(obj)
}
