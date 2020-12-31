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

#' Thanks to 
#' http://www.libpng.org/pub/png/book/chapter10.html
#' image_sample = light_out ^ gamma
#' it is said that
#' Once again, bear in mind that light_out and image_sample are scaled to the interval between 0 and 1;
#' that is, if the sample depth is 8 bits, the file samples range between 0 and 255, so image_sample is
#' obtained by dividing a given file sample by 255, in floating-point arithmetic. 
#' So,
#' image_sample = ymid and its range is [0,255]
#' light_out = xmid and its range is [xmin,xmax]
#' we have ymid / 255 = ((xmid - xmin) / (xmax - xmin)) ^ gamma
#' log(ymid / 255) = log((xmid - xmin) / (xmax - xmin)) * gamma
#' gamma = log(ymid / 255) / log((xmid - xmin) / (xmax - xmin))

#' @title Image Gamma Computation
#' @description
#' Computes image gamma transformation value.
#' @param V channel display properties containing 'xmin', 'xmax', 'xmid' and 'ymid'.
#' @keywords internal
computeGamma <- function(V) {
  if(missing(V)) stop("'V' can't be missing")
  stopifnot(all(c("xmin", "xmax", "xmid", "ymid") %in% names(V)))
  cpp_computeGamma(unlist(V[c("xmin", "xmax", "xmid", "ymid")]))
}

#' @title Image Gamma Modification
#' @description
#' Determines best xmid, ymid computes ymid for a given gamma
#' @param V named NumericVector of channel display properties containing 'xmin', 'xmax', 'xmid' and 'ymid'.
#' @param gamma gamma
#' @keywords internal
modifyGamma <- function(V, gamma = 1) {
  if(missing(V)) stop("'V' can't be missing")
  stopifnot(all(c("xmin", "xmax", "xmid", "ymid") %in% names(V)))
  O = unlist(V[c("xmin", "xmax", "xmid", "ymid")])
  f = function(O, gamma) { O["ymid"]= as.integer(255 * ((O["xmid"] - O["xmin"]) / (O["xmax"] - O["xmin"])) ^ gamma) }
  X = seq(O["xmin"], O["xmax"], 1)
  Y = sapply(X, FUN= function(x) {
    O["xmid"] <- x
    f(O, gamma)
  })
  fit = sapply(1:length(X), FUN=function(i) {
    O["xmid"] <- X[i]
    O["ymid"] <- Y[i]
    computeGamma(O)
  })
  best = which.min(fit - gamma)
  V["xmid"] <- X[best]
  V["ymid"] <- Y[best]
  V["gamma"] <- fit[best]
  return(V)
}
