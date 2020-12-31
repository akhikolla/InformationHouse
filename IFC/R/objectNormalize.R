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

#' @title Object Intensity Normalization
#' @description
#' Normalizes a matrix to [0,1].
#' @param mat a finite numeric matrix.
#' @param input_range a finite numeric vector of 2 values, sets the range of the input intensity values.\cr
#' Values exceeding this range are clipped. Default is c(0, 4095).
#' @param full_range if 'full_range' is TRUE, then 'input_range' will be set to c(0, 4095) and 'gamma' forced to 1. Default is FALSE.
#' @param force_range if 'force_range' is TRUE, then 'input_range' will be adjusted to object range in [-4095, +inf] and 'gamma' forced to 1. Default is FALSE.\cr
#' Note that this parameter takes the precedence over 'input_range' and 'full_range'.
#' @param gamma gamma correction. Default is 1, for no correction.
#' @details Note that negative values are used internally for removal of unmasked objects.
#' @return a [0,1] normalized matrix
#' @keywords internal
objectNormalize <- function(mat, input_range=c(0,4095), full_range=FALSE, force_range=FALSE, gamma=1) {
  cpp_normalize(mat = mat, input_range = input_range, full_range = full_range, force_range = force_range, gamma = gamma)
}
