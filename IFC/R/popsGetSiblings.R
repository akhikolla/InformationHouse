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

#' @title IFC_pops Sibling Population from Same Base Identification
#' @description
#' Gives names of graphical pops's siblings in a `IFC_data` object.
#' Siblings are built from a same population but different regions drawn on axes
#' of same feature(s) with same transformation(s) applied if any.
#' @param obj an `IFC_data` object extracted with features extracted.
#' @param pops graphical populations names to get siblings of.
#' @return names of population siblings.
#' @keywords internal
popsGetSiblings1 <- function(obj, pops) {
  if(missing(obj)) stop("'obj' can't be missing")
  if(missing(pops)) stop("'pops' can't be missing")
  if(!("IFC_data"%in%class(obj))) stop("'obj' is not of class `IFC_data`")
  if(length(obj$pops)==0) stop("please use argument 'extract_features' = TRUE with ExtractFromDAF() or ExtractFromXIF() and ensure that features were correctly extracted")
  if(is.null(pops)) stop("'pops' argument can't be NULL")
  N = names(obj$pops)
  if(!all(pops%in%N)) stop(paste0("pops [",paste0(pops[!pops%in%N], collapse = ","),"] not found in 'obj', valid names are:\n\t-", paste0(N, collapse="\n\t-")))
  lapply(obj$pops[pops], FUN=function(p) {
    if(p$type!="G") return(p$name)
    r = obj$regions[[p$region]]
    if(is.null(r)) return(p$name)
    map = sapply(obj$pops, FUN=function(m) {
      if(m$type!="G") return(rep(FALSE, 5))
      if(is.null(obj$regions[[m$region]])) return(rep(FALSE, 5))
      m$xlogrange = obj$regions[[m$region]]$xlogrange
      if(is.null(p$fy)) {
        return(c(m$base==p$base,
                 ifelse(is.null(m$fx), FALSE, identical(m$fx, p$fx)),
                 ifelse(is.null(m$xlogrange), FALSE, identical(r$xlogrange, m$xlogrange)), 
                 is.null(m$fy),
                 TRUE))
      } else {
        m$ylogrange = obj$regions[[m$region]]$ylogrange
        return(c(m$base==p$base,
                 ifelse(is.null(m$fx), FALSE, identical(m$fx, p$fx)), 
                 ifelse(is.null(m$xlogrange), FALSE, identical(m$xlogrange, r$xlogrange)), 
                 ifelse(is.null(m$fy), FALSE, identical(m$fy, p$fy)),
                 ifelse(is.null(m$ylogrange), FALSE, identical(m$ylogrange, r$ylogrange))))
      }
    })
    return(N[apply(map, 2, all)])
  })
}

#' @title IFC_pops Sibling Population from Same Region Identification
#' @description
#' Gives names of graphical pops's siblings in a `IFC_data` object.
#' Siblings are built from different populations but a same region drawn on axes
#' of same feature(s) with same transformation(s) applied if any.
#' @param obj an `IFC_data` object extracted with features extracted.
#' @param pops graphical populations names to get siblings of.
#' @return names of population siblings.
#' @keywords internal
popsGetSiblings2 <- function(obj, pops) {
  if(missing(obj)) stop("'obj' can't be missing")
  if(missing(pops)) stop("'pops' can't be missing")
  if(!("IFC_data"%in%class(obj))) stop("'obj' is not of class `IFC_data`")
  if(length(obj$pops)==0) stop("please use argument 'extract_features' = TRUE with ExtractFromDAF() or ExtractFromXIF() and ensure that features were correctly extracted")
  if(is.null(pops)) stop("'pops' argument can't be NULL")
  N = names(obj$pops)
  if(!all(pops%in%N)) stop(paste0("pops [",paste0(pops[!pops%in%N], collapse = ","),"] not found in 'obj', valid names are:\n\t-", paste0(N, collapse="\n\t-")))
  lapply(obj$pops[pops], FUN=function(p) {
    if(p$type!="G") return(p$name)
    r = obj$regions[[p$region]]
    if(is.null(r)) return(p$name)
    map = sapply(obj$pops, FUN=function(m) {
      if(m$type!="G") return(rep(FALSE, 5))
      if(is.null(obj$regions[[m$region]])) return(rep(FALSE, 5))
      m$xlogrange = obj$regions[[m$region]]$xlogrange
      if(is.null(p$fy)) {
        return(c(m$base!=p$base,
                 ifelse(is.null(m$fx), FALSE, identical(m$fx, p$fx)),
                 ifelse(is.null(m$xlogrange), FALSE, identical(r$xlogrange, m$xlogrange)), 
                 is.null(m$fy),
                 TRUE))
      } else {
        m$ylogrange = obj$regions[[m$region]]$ylogrange
        return(c(m$base!=p$base,
                 ifelse(is.null(m$fx), FALSE, identical(m$fx, p$fx)), 
                 ifelse(is.null(m$xlogrange), FALSE, identical(m$xlogrange, r$xlogrange)), 
                 ifelse(is.null(m$fy), FALSE, identical(m$fy, p$fy)),
                 ifelse(is.null(m$ylogrange), FALSE, identical(m$ylogrange, r$ylogrange))))
      }
    })
    return(N[apply(map, 2, all)])
  })
}

#' @title IFC_pops Sibling Population Identification
#' @description
#' Gives names of graphical pops's siblings in a `IFC_data` object.
#' @param obj an `IFC_data` object extracted with features extracted.
#' @param pops graphical populations names to get siblings of.
#' @return names of population siblings.
#' @keywords internal
popsGetSiblings <- function(obj, pops) {
  foo = c()
  bar = unique(c(unlist(popsGetSiblings1(obj = obj, pops = pops)), unlist(popsGetSiblings2(obj = obj, pops = pops))))
  while(!identical(foo, bar)) {
    foo = bar
    bar = unique(c(unlist(popsGetSiblings1(obj = obj, pops = bar)), unlist(popsGetSiblings2(obj = obj, pops = bar))))
  }
  type = sapply(obj$pops[bar], FUN = function(p) p$type)
  if(!all(type == "G")) return(NULL)
  x = sapply(obj$pops[bar], FUN = function(p) p$fx)
  type = sapply(obj$pops[bar], FUN = function(p) obj$regions[[p$region]]$type == "line")
  if(length(x) != length(type)) return(NULL)
  if((length(unique(x)) != 1) ||
     (length(unique(type)) != 1)) return(NULL)
  xlog = sapply(obj$pops[bar], FUN = function(p) obj$regions[[p$region]]$xlogrange)
  if(!all(type)) {
    y = unlist(lapply(obj$pops[bar], FUN = function(p) p$fy))
    ylog = sapply(obj$pops[bar], FUN = function(p) obj$regions[[p$region]]$ylogrange)
    if((length(unique(y)) != 1) ||
       (length(unique(xlog)) != 1) ||
       (length(unique(ylog)) != 1)) return(NULL)
  } else {
    if(length(unique(xlog)) != 1) return(NULL)
  }
  return(bar)
}