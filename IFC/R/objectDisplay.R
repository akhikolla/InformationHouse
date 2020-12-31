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

#' @title Object Display
#' @description
#' This function is intended to display object extracted by \code{\link{objectExtract}}.
#' @param image An object extracted by \code{\link{objectExtract}} of class `IFC_img` or `IFC_msk`.\cr
#' Note that a matrix with finite values can also be used.
#' @param input_range a finite numeric vector of 2 values, sets the range of the input intensity values.\cr
#' Values exceeding this range are clipped. Default is 'c(0, 4095)'.
#' @param full_range if 'full_range' is TRUE, then 'input_range' will be set to 'c(0, 4095)' and 'gamma' forced to 1. Default is FALSE.
#' @param force_range if 'force_range' is TRUE, then 'input_range' will be adjusted to object range in [-4095, +inf] and 'gamma' forced to 1. Default is FALSE.\cr
#' Note that this parameter takes the precedence over 'input_range' and 'full_range'.
#' @param gamma gamma correction. Default is 1, for no correction.
#' @param color a color. Default is "Green".
#' @details If input 'image' is of class `IFC_img` or `IFC_msk`, then if 'input_range', 'full_range', 'force_range', 'gamma' and / or 'color' parameters is/are missing,
#' it/they will be extracted from 'image' attributes.\cr
#' If input 'image' is not of one of class `IFC_img` or `IFC_msk`, then force_range will be forced to TRUE.\cr
#' An error will be thrown if input image contains non finite values.
#' @param dpi display resolution. Default is 300.
#' @return it invisibly returns NULL
#' @export
objectDisplay = function(image, input_range = c(0, 4095), full_range = FALSE, force_range = FALSE, gamma = 1, color = "Green", dpi = 300) {
  dpi = na.omit(as.integer(dpi)); dpi = dpi[dpi>0]; dpi = dpi[is.finite(dpi)]
  assert(dpi, len = 1, typ = "integer")
  d = dim(image)
  K = class(image)
  foo = NULL
  if(any(c("IFC_img", "IFC_msk") %in% K)) {
    if(missing(input_range)) input_range = attr(image, "input_range")
    if(missing(full_range)) force_range = attr(image, "full_range")
    if(missing(force_range)) force_range = attr(image, "force_range")
    if(missing(gamma)) gamma = attr(image, "gamma")
    if(missing(color)) color = attr(image, "color")
    checkColor(color)
    if("IFC_img" %in% K) {
      foo = objectColorize(objectNormalize(attr(image, "raw"), input_range = input_range, full_range = full_range, force_range = force_range, gamma = gamma), color)
    } else {
      foo = objectColorize(objectNormalize(attr(image, "raw"), input_range = input_range, full_range = full_range, force_range = force_range, gamma = 1), color)
    }
  } else {
    if("matrix" %in% K) {
      checkColor(color)
      foo = objectColorize(objectNormalize(image, force_range = TRUE), color)
    }
  }
  if(length(foo) == 0) stop("'image' is not compatible with objectDisplay")
  grid.newpage()
  do.call(what = "grid.raster", args = list(image = foo,
                                            width = unit(dpi * d[2] / 96, "points"),
                                            height = unit(dpi * d[1] / 96, "points"),
                                            interpolate = FALSE))
}
