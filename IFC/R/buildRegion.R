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

#' @title IFC Region Coercion
#' @description
#' Helper to build a list to allow region export.
#' @param type Region's type. Either "line", "rect", "poly" or "oval".
#' @param label label of the region.
#' @param cx x label's position. If not provided x center will be used.
#' @param cy y label's position. If not provided y center will be used.
#' @param color color of the region. See \code{\link{paletteIFC}} for allowed colors.
#' @param lightcolor lightcolor of the region. See \code{\link{paletteIFC}} for allowed colors.
#' @param ismarker Default is 'false'. Allowed are 'true' or 'false'. Used for compatibility with amnis file but role remains unknown.
#' @param doesnotoverride Default is 'false'. Allowed are 'true' or 'false'. Used for compatibility with amnis file but role remains unknown.
#' @param xlogrange determines hyper parameter of smoothLinLog transformation for x-axis. Default is "P" for no transformation.
#' @param ylogrange determines hyper parameter of smoothLinLog transformation for y-axis. Default is "P" for no transformation.
#' @param x vector of x vertices values.
#' @param y vector of y vertices values.
#' @param ... Other arguments to be passed.
#' @return a list containing all region information.
#' @export
buildRegion <- function(type, label, cx, cy, color, lightcolor, ismarker="false", doesnotoverride="false", xlogrange, ylogrange, x, y, ...) {
  dots = list(...)
  if(missing(label)) stop("'label' can't be missing")
  assert(label, len=1, typ="character")
  if(grepl("|",label,fixed=TRUE)) warning(paste0("'|' found in 'label': ",label, ", this may cause unstable behaviour"), call. = FALSE, immediate. = TRUE)
  if(missing(type)) stop("'type' can't be missing")
  assert(type, len=1, alw=c("line","rect","poly","oval"))
  ismarker = tolower(ismarker); assert(ismarker, len=1, alw=c("false","true"))
  doesnotoverride = tolower(doesnotoverride); assert(doesnotoverride, len=1, alw=c("false","true"))
  if(missing(xlogrange)) stop("'xlogrange' can't be missing")
  if(missing(ylogrange)) stop("'ylogrange' can't be missing")
  xlogrange = as.character(xlogrange)
  ylogrange = as.character(ylogrange)
  if(xlogrange != "P") {
    tmp = as.integer(xlogrange)
    if(is.na(tmp)) stop("'xlogrange' should be coercible to positive integer")
    if(tmp<0) stop("'xlogrange' should be coercible to positive integer")
    xlogrange = as.character(tmp)
  }
  if(ylogrange != "P") {
    tmp = as.integer(ylogrange)
    if(is.na(tmp)) stop("'ylogrange' should be coercible to positive integer")
    if(tmp<0) stop("'ylogrange' should be coercible to positive integer")
    ylogrange = as.character(tmp)
  }
  if(missing(color)) {
    if(missing(lightcolor)) {
      tmp = sample(nrow(paletteIFC("")),1)
    } else {
      assert(lightcolor, len=1, alw=unlist(paletteIFC("")))
      tmp = which(paletteIFC("")%in%lightcolor, arr.ind=TRUE)[1]
      if(is.na(tmp)) tmp = sample(nrow(paletteIFC("")),1) 
    }
    color = paletteIFC("")$color[tmp]
    lightcolor = paletteIFC("")$lightModeColor[tmp]
  } else {
    if(color%in%paletteIFC("")$color_R) color = paletteIFC("")$color[color==paletteIFC("")$color_R][1]
    assert(color, len=1, alw=paletteIFC("palette"))
  }
  if(missing(lightcolor)) {
    if(missing(color)) {
      tmp = sample(nrow(paletteIFC("")),1)
    } else {
      assert(color, len=1, alw=unlist(paletteIFC("")))
      tmp = which(color==paletteIFC(""), arr.ind=TRUE)[1]
      if(is.na(tmp)) tmp = sample(nrow(paletteIFC("")),1)
    }
    color = paletteIFC("")$color[tmp]
    lightcolor = paletteIFC("")$lightModeColor[tmp]
  } else {
    if(lightcolor%in%paletteIFC("")$lightModeColor_R) lightcolor = paletteIFC("")$lightModeColor[lightcolor==paletteIFC("")$lightModeColor_R][1]
    assert(lightcolor, len=1, alw=paletteIFC("palette"))
  }
  
  if(missing(x)) stop("'x' can't be missing")
  if(missing(y)) stop("'y' can't be missing")
  x = as.numeric(x)
  y = as.numeric(y)
  x = x[is.finite(x)]
  y = y[is.finite(y)]
  if(type=="poly") {
    if(length(x) != length(y)) stop("'x' and 'y' should be numeric vectors of equal length with finite values")
    if(length(x)<2) stop("type='poly' and number of vertices is smaller than 2")
    if(missing(cx) | missing(cy)) {
      # cent = cpp_poly_centroid(cbind(x,y)) # TODO # for future use
      # if(missing(cx)) cx= cent[1]
      # if(missing(cy)) cy= cent[2]
      if(missing(cx)) cx= mean(x)
      if(missing(cy)) cy= mean(y)
    }
  } else {
    if(length(x) != 2) stop(paste0("'x' should be a length 2 numeric vector with finite values when type is '",type,"'"))
    if(length(y) != 2) stop(paste0("'y' should be a length 2 numeric vector with finite values when type is '",type,"'"))
    if(type=="line" & y[1]!=y[2]) stop("'y' values should be identical when type is 'line")
    if(missing(cx) | missing(cy)) {
      if(missing(cx)) cx= mean(x)
      if(missing(cy)) cy= mean(y)
    }
  }
  return(list("type"=type, "label"=label, "cx"=cx, "cy"=cy, "color"=color, "lightcolor"=lightcolor, "ismarker"=ismarker, "doesnotoverride"=doesnotoverride, "xlogrange"=xlogrange, "ylogrange"=ylogrange, "x"=x, "y"=y))
}
